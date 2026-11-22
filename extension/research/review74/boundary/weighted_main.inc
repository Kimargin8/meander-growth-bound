// Find and certify a positive vector for the restricted exit language.
// Arguments: original stage-73 checkpoint, additional sweeps, new checkpoint.
int main(int argc, char** argv) {
    if (argc != 4 || std::string(argv[1]) == argv[3]) return 2;
    std::ifstream input(argv[1], std::ios::binary);
    std::array<U, 5> header{};
    input.read(reinterpret_cast<char*>(header.data()), sizeof(header));
    if (!input || header[0] != 0x5354524950373301ULL || header[1] < 1 || header[1] > 16)
        throw std::runtime_error("invalid starting checkpoint");
    const int extra = std::stoi(argv[2]);
    if (extra < 0 || extra > 2000) return 2;
    Board board(2 * header[1]);
    int round = 0, shift = 0;
    board.checkpoint(argv[1], true, round, shift);
    const int end = round + extra;
    for (; round <= end; ++round) {
        if (round == end) board.checkpoint(argv[3], false, round, shift);
        board.sweep<false, false>(0, shift, 0);
        bool measure = (end - round) % 5 == 0 || round == int(header[3]);
        Stats result = measure ? board.sweep<false, true>(1, 0, shift)
                               : board.sweep<false, false>(1, 0, shift);
        if (measure)
            std::cout << "round " << round << " min " << result.lo.a << "/" << result.lo.b
                      << " max " << result.hi.a << "/" << result.hi.b << std::endl;
        shift = 0;
        while ((result.top >> shift) > scale) ++shift;
    }
}
