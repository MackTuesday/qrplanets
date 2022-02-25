#include <random>
#include "BrentsRand.h"

static std::random_device rd;

static unsigned lastRand = rd();

void BrentsRandSet(unsigned seed) {
    lastRand = seed;
}
unsigned BrentsRand(unsigned max) {
    lastRand = ((lastRand * 1103515245) + 12345) % 2147483647;
    return ((lastRand >> 15) * max) >> 16;
}
double BrentsUnitRand() {
    lastRand = ((lastRand * 1103515245) + 12345) % 2147483647;
    return (lastRand >> 15) / 65536.0;
}
