//
// Created by benjamin on 2/11/22.
//

#include <cstdio>
#include <string>
#include "def.h"


int main(int argc, char **argv) {
    if (argc == 1) {
        printf("Need to give argument");
        return 0;
    }

    int test_attr = atoi(argv[1]);

    std::string keyInput = "";
    std::string encInput = "";

    uint32_t N_ATTR = test_attr;

    uint32_t *attr_int_list = NULL;
    attr_int_list = (uint32_t *) malloc(sizeof(uint32_t) * test_attr);

    int d = 1;

    for (int k = 0; k < test_attr; k++) {
        keyInput = keyInput + "attr" + std::to_string(d);
        encInput = encInput + "attr" + std::to_string(d);

        if (k < test_attr - 1) {
            keyInput = keyInput + "|";
            encInput = encInput + " and ";
        }

        attr_int_list[k] = d;

        d++;

    }

    struct master_key msk;
    init_master_key(N_ATTR, &msk);

    return 0;
}
