#include "lib/policy/policy_tree.h"
int main(int argc, char const *argv[])
{
    test1();
    return 0;
}

int test1(){
    char formula[] = "AND(OR(A), OR(B), OR(OR(C), OR(D))";
    tree_from_string(formula);
    
}