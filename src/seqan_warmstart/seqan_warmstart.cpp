// seqan_warmstart.cpp
//
// Heuristic-MSA "warm start" aligner: a drop-in replacement for the external
// MUSCLE invocation in run_MUSCLE(). It accepts MUSCLE-compatible command-line
// flags, so the calling command string only needs its binary name changed. It
// reads a FASTA, runs a SeqAn v2 progressive multiple sequence alignment under
// a user-supplied substitution matrix and affine gap penalties, and writes an
// aligned FASTA.
//
// The downstream solver re-scores the alignment with its own sum-of-pairs
// objective, so this tool only needs to emit a *valid* alignment whose
// sequence headers match the input verbatim. Alignment quality affects only
// the strength of the initial incumbent (hence solver speed), never the
// correctness of the exact algorithm.
//
// SeqAn v2 (2.5.x); built at C++17. See the `seqan_warmstart` Makefile target.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/score.h>
#include <seqan/seq_io.h>

using namespace seqan2;

namespace {

typedef String<char>                            TSequence;
typedef Score<int, ScoreMatrix<char, Default> > TScore;

// All substitution-matrix symbols and sequence residues are 7-bit ASCII, so we
// only ever setScore() over this range. SeqAn's setScore casts the symbol to
// (signed) char and then to unsigned int; a byte >= 128 therefore becomes a
// negative char and a huge index that corrupts memory on platforms where char
// is signed (x86) -- harmless where char is unsigned (ARM), which masks the bug
// in testing. Restricting to ASCII keeps every index non-negative and portable.
const unsigned ASCII_SIZE = 128;

// Parsed command line. Defaults mirror the flag semantics MUSCLE used under
// the caller (a length-L gap costs gapOpen + (L-1)*gapExtend in both MUSCLE and
// SeqAn, so the flag values map across directly).
struct Options {
    std::string inPath;
    std::string outPath;
    std::string matrixPath;
    int  gapOpen   = -1;   // score of the first character of a gap
    int  gapExtend = -1;   // score of each subsequent gap character
    bool quiet     = false;
};

// Round a possibly-fractional CLI value (e.g. "-14.220000") to int. SeqAn's
// Score<int, ...> is integral. For the morphology use case the penalties are
// already integers (lossless); only the protein defaults round (e.g. -14.22 ->
// -14). See session notes for the empirical-consequence discussion.
int roundFlag(const std::string & s) {
    return static_cast<int>(std::lround(std::atof(s.c_str())));
}

// Parse MUSCLE-style flags. Recognised value flags consume the next token;
// recognised boolean flags consume nothing. Unsupported flags (-center and the
// -sp family, etc.) are ignored silently as required -- -center additionally
// consumes its value so it is not mistaken for a positional token.
bool parseArgs(int argc, char ** argv, Options & opt) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto value = [&](const char * name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "seqan_warmstart: missing value for " << name << "\n";
                return std::string();
            }
            return std::string(argv[++i]);
        };
        if      (a == "-in")        opt.inPath     = value("-in");
        else if (a == "-out")       opt.outPath    = value("-out");
        else if (a == "-matrix")    opt.matrixPath = value("-matrix");
        else if (a == "-gapopen")   opt.gapOpen    = roundFlag(value("-gapopen"));
        else if (a == "-gapextend") opt.gapExtend  = roundFlag(value("-gapextend"));
        else if (a == "-quiet")     opt.quiet      = true;
        else if (a == "-refine")    { /* SeqAn has no in-place refine: re-align */ }
        else if (a == "-center")    { (void) value("-center"); } // ignored value flag
        else if (!a.empty() && a[0] == '-') { /* unsupported flag: ignore */ }
        else { /* stray positional token: ignore */ }
    }
    if (opt.inPath.empty() || opt.outPath.empty()) {
        std::cerr << "seqan_warmstart: both -in and -out are required\n";
        return false;
    }
    return true;
}

// Load an NCBI-format substitution matrix, parsed exactly as the downstream
// solver's own loader does so the two stay in agreement:
//   * first non-comment row: whitespace-separated column header characters;
//   * each subsequent row: a row-label character then integer scores.
// Out-of-alphabet pairs are filled with the most-negative declared score (the
// ScoreMatrix Default otherwise zero-fills, which could beat a real mismatch).
// Valid input never references an undeclared character, so this is defensive.
bool loadMatrix(const std::string & path, TScore & sc) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        std::cerr << "seqan_warmstart: cannot open matrix file: " << path << "\n";
        return false;
    }
    std::vector<char> cols;
    std::vector<std::pair<char, std::vector<int> > > rows;
    int  minScore  = 0;
    bool haveScore = false;
    bool header    = true;
    std::string line;
    while (std::getline(in, line)) {
        bool blank = true;
        for (std::string::size_type k = 0; k < line.size(); ++k)
            if (!std::isspace((unsigned char) line[k])) { blank = false; break; }
        if (blank || line[0] == '#') continue;       // skip blanks and comments
        std::istringstream ss(line);
        std::string tok;
        if (header) {
            while (ss >> tok) cols.push_back(tok[0]);
            header = false;
        } else {
            if (!(ss >> tok)) continue;
            char rowCh = tok[0];
            std::vector<int> scores;
            int v;
            while (ss >> v) {
                scores.push_back(v);
                if (!haveScore || v < minScore) { minScore = v; haveScore = true; }
            }
            rows.push_back(std::make_pair(rowCh, scores));
        }
    }
    if (cols.empty() || rows.empty()) {
        std::cerr << "seqan_warmstart: empty or invalid matrix file: " << path << "\n";
        return false;
    }
    for (unsigned i = 0; i < ASCII_SIZE; ++i)
        for (unsigned j = 0; j < ASCII_SIZE; ++j)
            setScore(sc, (char) i, (char) j, minScore);
    for (std::size_t r = 0; r < rows.size(); ++r)
        for (std::size_t c = 0; c < rows[r].second.size() && c < cols.size(); ++c)
            setScore(sc, rows[r].first, cols[c], rows[r].second[c]);
    return true;
}

} // namespace

int main(int argc, char ** argv) {
    Options opt;
    if (!parseArgs(argc, argv, opt)) return 1;

    // --- Read input FASTA (SeqAn native). ids retain the full header line,
    //     including any spaces, which the solver's name lookup depends on.
    StringSet<CharString> ids;
    StringSet<TSequence>  seqs;
    {
        // Open via an explicit stream so SeqAn detects the format from file
        // *content* (the leading '>'), not the filename extension. The default
        // path-based open() rejects anything not named .fa/.fasta/...; MSA-MDD
        // only ever passes ".fa", but this keeps the tool robust to any -in
        // path (e.g. a "_s.txt" alignment fed back for refinement).
        std::ifstream ifs(opt.inPath.c_str(), std::ios::binary);
        if (!ifs.good()) {
            std::cerr << "seqan_warmstart: cannot open input: " << opt.inPath << "\n";
            return 1;
        }
        SeqFileIn fin;
        if (!open(fin, ifs)) {
            std::cerr << "seqan_warmstart: unrecognised or empty FASTA input: "
                      << opt.inPath << "\n";
            return 1;
        }
        try {
            readRecords(ids, seqs, fin);
        } catch (Exception const & e) {
            std::cerr << "seqan_warmstart: FASTA read error: " << e.what() << "\n";
            return 1;
        }
    }

    // Strip gap characters from input sequences. The refine call path feeds an
    // already-aligned (gapped) FASTA; SeqAn would otherwise treat '-' as a
    // residue. Stripping makes both the fresh and refine paths re-align from
    // ungapped sequences (SeqAn exposes no in-place refinement).
    for (unsigned i = 0; i < length(seqs); ++i) {
        TSequence cleaned;
        for (unsigned j = 0; j < length(seqs[i]); ++j)
            if (seqs[i][j] != '-') appendValue(cleaned, seqs[i][j]);
        seqs[i] = cleaned;
    }

    // --- Build the scoring scheme.
    TScore sc;
    if (!opt.matrixPath.empty()) {
        if (!loadMatrix(opt.matrixPath, sc)) return 1;
    } else {
        // No matrix supplied: simple identity scheme (match +1, mismatch -1).
        for (unsigned i = 0; i < ASCII_SIZE; ++i)
            for (unsigned j = 0; j < ASCII_SIZE; ++j)
                setScore(sc, (char) i, (char) j, (i == j) ? 1 : -1);
    }
    setScoreGapOpen(sc,   opt.gapOpen);    // first gap char  (MUSCLE -gapopen)
    setScoreGapExtend(sc, opt.gapExtend);  // later gap chars (MUSCLE -gapextend)

    // --- Align. globalMsaAlignment needs >= 2 sequences; pass trivial inputs
    //     (0 or 1 sequence) straight through.
    Align<TSequence> aln;
    resize(rows(aln), length(seqs));
    for (unsigned i = 0; i < length(seqs); ++i)
        assignSource(row(aln, i), seqs[i]);
    bool aligned = (length(seqs) >= 2);
    if (aligned)
        globalMsaAlignment(aln, sc);

    // --- Write aligned FASTA. Reproduce each input header verbatim (prefixed
    //     with '>') so the solver's Load_sol maps rows back to sequences.
    std::ofstream out(opt.outPath.c_str());
    if (!out.good()) {
        std::cerr << "seqan_warmstart: cannot open output: " << opt.outPath << "\n";
        return 1;
    }
    typedef Row<Align<TSequence> >::Type TRow;
    for (unsigned i = 0; i < length(seqs); ++i) {
        out << '>' << ids[i] << '\n';
        if (aligned) {
            TRow & r = row(aln, i);
            for (unsigned col = 0; col < length(r); ++col)
                out << (isGap(r, col) ? '-' : (char) value(r, col));
        } else {
            out << seqs[i];
        }
        out << '\n';
    }
    out.close();

    if (!opt.quiet)
        std::cerr << "seqan_warmstart: wrote " << length(seqs)
                  << " sequences to " << opt.outPath << "\n";
    return 0;
}
