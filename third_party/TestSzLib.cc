#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C"{
#endif
#include "szlib.h"
#ifdef __cplusplus
}
#endif

#define OPTIONS_MASK (SZ_RAW_OPTION_MASK \
                      | SZ_MSB_OPTION_MASK \
                      | SZ_NN_OPTION_MASK)
#define PIXELS_PER_BLOCK (8)
#define PIXELS_PER_SCANLINE (PIXELS_PER_BLOCK*128)

void createFile(const char* fileName, const char* fileContent){
    FILE *fptr;

    // Open a file in writing mode
    fptr = fopen(fileName, "w");

    // Write some text to the file
    fprintf(fptr, fileContent);

    // Close the file
    fclose(fptr);
}

int main(int argc, char *argv[])
{
    printf("Write temporary file...");

    const char* firstFileName = "/tmp/test_sz_lib_temp_file.txt";
    const char* secondFileName = "/tmp/test_sz_lib_temp_file2.txt";

    createFile(firstFileName,"test sz compress random 1234567890-.");
    createFile(secondFileName,"test another sz compress.");

    printf("Starting...");
    int status = 0;
    SZ_com_t sz_param;
    unsigned char *source = NULL;
    unsigned char *dest = NULL;
    unsigned char *dest1 = NULL;
    size_t destLen, dest1Len, sourceLen;
    FILE *fp;

    if ((fp = fopen(firstFileName, "rb")) == NULL) {
        fprintf(stderr, "Can't open %s\n", argv[1]);
        return 1;
    }

    fseek(fp, 0L, SEEK_END);
    sourceLen = ftell(fp);
    fseek(fp, 0L, SEEK_SET);
    destLen = sourceLen + sourceLen / 10;

    sz_param.options_mask = OPTIONS_MASK;
    sz_param.bits_per_pixel = 64;
    sz_param.pixels_per_block = PIXELS_PER_BLOCK;
    sz_param.pixels_per_scanline = PIXELS_PER_SCANLINE;

    source = (unsigned char *)malloc(sourceLen);
    dest = (unsigned char *)malloc(destLen);
    dest1 = (unsigned char *)malloc(destLen);

    if (source == NULL || dest == NULL || dest1 == NULL) {
        status = 99;
        goto DESTRUCT;
    }

    sourceLen = fread(source, 1, sourceLen, fp);

    status = SZ_BufftoBuffCompress(dest, &destLen,
                                   source, sourceLen, &sz_param);
    if (status != SZ_OK)
        goto DESTRUCT;

    dest1Len = sourceLen;
    status = SZ_BufftoBuffDecompress(dest1, &dest1Len,
                                     dest, destLen, &sz_param);
    if (status != SZ_OK)
        goto DESTRUCT;

    if (memcmp(source, dest1, sourceLen) != 0)
        fprintf(stderr, "File %s Buffers differ\n", argv[2]);

DESTRUCT:
    if (source)
        free(source);
    if (dest)
        free(dest);
    if (dest1)
        free(dest1);
    printf("Status = %d", status);
    return status;
}