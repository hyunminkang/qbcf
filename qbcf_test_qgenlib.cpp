#include <vector>
#include "qgenlib/params.h"
#include "qgenlib/qgen_error.h"
int qbcf_test_qgenlib(int argc, char** argv) {
    std::string str1;
    std::vector<std::string> str2;
    int int3 = 0;
    double double4 = 0.0;
    bool bool5 = false;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
        LONG_PARAM_GROUP("Options for string arguments", NULL)
        LONG_STRING_PARAM("str1", &str1, "A dummy string argument")
        LONG_MULTI_STRING_PARAM("str2", &str2, "Another dummy string argument (can be specified multiple times)")

        LONG_PARAM_GROUP("Options for non-string arguments", NULL)
        LONG_INT_PARAM("int3", &int3, "A dummy integer argument")    
        LONG_DOUBLE_PARAM("double4", &double4, "A dummy double argument")
        LONG_PARAM("bool5",&bool5,"A dummy boolean argument")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    printf("str1: %s\n", str1.c_str());
    printf("str2: ");
    for (int i = 0; i < str2.size(); i++) {
        printf("%s ", str2[i].c_str());
    }
    printf("\n");
    printf("int3: %d\n", int3);
    printf("double4: %f\n", double4);
    printf("bool5: %s\n", bool5 ? "true" : "false");

    notice("This is an example notice message");
    warning("This is an example warning message");

    return 0;
}