#include "util.h"
#include <cmath>
#include <map>

#define MIN(x, y) (x < y ? x : y)

num_t distance(vec_t& a, vec_t& b) {
    num_t x = a.x - b.x;
    num_t y = a.y - b.y;
    return std::round(std::sqrt(x * x + y * y));
}

void compute_neighbours(int n, city_t *world, city_t& city) {
    // lexicographic sorting on the key
    std::map<std::pair<num_t, int>, city_t*> collection;
    for (int i = 0; i < n; ++i) {
        city_t* other = &world[i];
        // don't add the city itself
        if (other == &city) {
            continue;
        }

        num_t dist = distance(other->pos, city.pos);
        // we add the index as a second value, otherwise two equidistant cities would clash (omitting one)
        collection.insert(std::pair(std::pair(dist, i), other));
        // FIXME use another implementation as this currently has a nasty runtime
    }

    int i = 0;
    for(const auto & [key, value] : collection) {
        if (i >= NEAREST) {
            break;
        }
        city.nearest[i] = value;
        i++;
    }
}
