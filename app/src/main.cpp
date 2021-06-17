#define DOCTEST_CONFIG_IMPLEMENT

#include <doctest/doctest.h>
#include <iostream>

int main(int argc, char** argv) {
    doctest::Context context;
    context.applyCommandLine(argc, argv);

    //int res = context.run();

    // if (context.shouldExit()) {
    //     return res;
    // }

    std::cout << "Hello cherif" << '\n';

    return 1;
}
