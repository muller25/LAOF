#include "Image.h"

#include <vector>

int main(int argc, char *argv[])
{
    DImage src(10, 10);
    DImage dst(src);
    DImage eq = src;
    std::vector<DImage> v;

    printf("start to push...");
    v.push_back(src);
    v.push_back(dst);
    v.push_back(eq);
    printf("done\n");
    
    return 0;
}
