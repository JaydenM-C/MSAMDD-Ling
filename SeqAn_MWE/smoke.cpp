// SeqAn v2 standalone smoke test (Phase 1).
//
// Purpose: prove that the vendored SeqAn 2.5.3 tree can drive an exact-API
// progressive MSA (globalMsaAlignment) over an ARBITRARY ASCII alphabet using a
// user-supplied substitution matrix -- the capability MUSCLE lacks and the whole
// reason SeqAn is replacing it as MSA-MDD's warm-start aligner.
//
// Note: SeqAn is vendored as the full upstream repo (git submodule), so headers
// live at includes/seqan/include/seqan/... and the include root is
// includes/seqan/include -- NOT includes/. See build.sh.
//
// Build (from this folder):  ./build.sh
//   equivalently: g++ -std=gnu++17 -I ../includes/seqan/include smoke.cpp -o smoke
// Run:                       ./smoke tiny.fa

#include <iostream>
#include <seqan/seq_io.h>      // SeqFileIn, readRecords
#include <seqan/align.h>       // Align, rows, assignSource
#include <seqan/graph_msa.h>   // globalMsaAlignment

using namespace seqan2;

int main(int argc, char ** argv)
{
    const char * path = (argc > 1) ? argv[1] : "tiny.fa";

    // --- Read a tiny FASTA over a raw char alphabet (no DNA/protein assumption).
    StringSet<CharString>     ids;
    StringSet<String<char> >  seqs;
    SeqFileIn fin;
    if (!open(fin, path)) { std::cerr << "ERROR: cannot open " << path << "\n"; return 1; }
    readRecords(ids, seqs, fin);
    std::cout << "Read " << length(seqs) << " sequences from " << path << "\n\n";

    // --- Custom substitution matrix: positive diagonal, negative off-diagonal,
    //     affine gaps. (Mirrors the morphology scoring scheme in the project
    //     conventions: match rewarded, mismatch penalised, no natural classes.)
    typedef Score<int, ScoreMatrix<char, Default> > TScore;
    TScore sc;
    for (unsigned i = 0; i < ValueSize<char>::VALUE; ++i)
        for (unsigned j = 0; j < ValueSize<char>::VALUE; ++j)
            setScore(sc, (char)i, (char)j, (i == j) ? 2 : -1);
    setScoreGapOpen(sc,   -2);   // affine gap-open
    setScoreGapExtend(sc, -1);   // affine gap-extend

    // --- Run the exact SeqAn API used downstream by seqan_warmstart.
    Align<String<char> > aln;
    resize(rows(aln), length(seqs));
    for (unsigned i = 0; i < length(seqs); ++i)
        assignSource(row(aln, i), seqs[i]);

    globalMsaAlignment(aln, sc);

    std::cout << aln << std::endl;
    return 0;
}
