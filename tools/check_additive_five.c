/* Independent exhaustive finite filter, not a Euclidean chromatic certificate.
 * Input: text lines, each containing one 12-coordinate F5 direction vector.
 * Required: all twelve coordinate vectors occur. A proper F5-linear coloring
 * must then have every coefficient nonzero. Scale coefficient zero to 1 and
 * enumerate the remaining 4^11 choices. No SAT or hn imports are used.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: %s directions.txt\n", argv[0]); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("directions"); return 2; }
    unsigned directions[1024][12], count = 0, basis = 0;
    while (1) {
        unsigned row[12]; int read = fscanf(f, "%u", &row[0]);
        if (read == EOF) break;
        if (read != 1 || count == 1024) return 2;
        for (unsigned j=1; j<12; ++j) if (fscanf(f,"%u",&row[j]) != 1) return 2;
        unsigned support=0, position=0;
        for (unsigned j=0;j<12;++j) {
            if(row[j]>4) return 2;
            if(row[j]) {++support;position=j;}
        }
        if(!support) return 2;
        if(support==1) { basis|=1u<<position; continue; }
        for(unsigned j=0;j<12;++j) directions[count][j]=row[j];
        ++count;
    }
    fclose(f);
    if(basis != 4095) {fputs("missing coordinate directions\n",stderr);return 2;}
    uint64_t checked=0, survivors=0;
    for(uint32_t assignment=0; assignment<(1u<<22); ++assignment) {
        unsigned coefficients[12]={1}; uint32_t digits=assignment;
        for(unsigned j=1;j<12;++j) {coefficients[j]=1+(digits&3);digits>>=2;}
        unsigned d;
        for(d=0;d<count;++d) {
            unsigned sum=0;
            for(unsigned j=0;j<12;++j)sum+=directions[d][j]*coefficients[j];
            if(sum%5==0)break;
        }
        ++checked;
        if(d==count) {
            ++survivors;
            if(survivors==1) {
                fputs("first_survivor",stdout);
                for(unsigned j=0;j<12;++j)printf(" %u",coefficients[j]);
                putchar('\n');
            }
        }
    }
    printf("assignments_checked %llu\nsurvivors %llu\n",
           (unsigned long long)checked,(unsigned long long)survivors);
    puts("scope additive F5 functionals only; no Euclidean graph lower bound");
    return survivors ? 1 : 0;
}
