/***************************************************************************
 *  Description:
 *      Generate a matrix of allelic-depths from single-sample VCF files
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-09  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include <limits.h>
#include <tsvio.h>
#include "ad-matrix.h"

int     main(int argc,char *argv[])

{
    file_list_t file_list;
    char        *list_filename = argv[1],
		*matrix_filename = argv[2];
    
    if ( argc != 3 )
	usage(argv);
    open_files(list_filename, &file_list, "r");
    build_matrix(&file_list, matrix_filename);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Read a list of VCF files from filename
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-09  Jason Bacon Begin
 ***************************************************************************/

void    open_files(char *list_filename, file_list_t *file_list, char *mode)

{
    FILE        *fp;
    extern int  errno;
    char        *temp_filename;
    size_t      actual_len,
		c;
    
    if ( (fp = fopen(list_filename, "r")) == NULL )
    {
	fprintf(stderr, "ad-matrix: Cannot open %s: %s\n",
		list_filename, strerror(errno));
	exit(EX_DATAERR);
    }
    
    if ( (temp_filename = malloc(PATH_MAX + 1)) == NULL )
    {
	fprintf(stderr, "open_files(): Cannot allocate temp filename.\n");
	exit(EX_UNAVAILABLE);
    }
    
    // Count VCF filenames
    file_list->count = 0;
    while ( fgets(temp_filename, PATH_MAX, fp) != NULL )
	++file_list->count;
    printf("%zu VCF files.\n", file_list->count);
    
    // Allocate list and read VCF filenames
    file_list->filename = (char **)malloc(file_list->count * sizeof(char *));
    if ( file_list == NULL )
    {
	fprintf(stderr, "open_files(): Cannot allocate array.\n");
	exit(EX_UNAVAILABLE);
    }
    file_list->fp = (FILE **)malloc(file_list->count * sizeof(FILE *));
    if ( file_list == NULL )
    {
	fprintf(stderr, "open_files(): Cannot allocate array.\n");
	exit(EX_UNAVAILABLE);
    }
    rewind(fp);
    for (c = 0; c < file_list->count; ++c)
    {
	tsv_read_field(fp, temp_filename, PATH_MAX, &actual_len);
	if ( (file_list->filename[c] = strdup(temp_filename)) == NULL )
	{
	    fprintf(stderr,
		    "open_files(): Error allocating filename[%zu]\n", c);
	    exit(EX_UNAVAILABLE);
	}
	if ( (file_list->fp[c] = fopen(file_list->filename[c], mode)) == NULL )
	{
	    fprintf(stderr, "open_file_list(): Cannot open %s: %s\n",
		    file_list->filename[c], strerror(errno));
	    exit(EX_UNAVAILABLE);   // FIXME: Tailor to mode?
	}
    }        
    fclose(fp);
}


/***************************************************************************
 *  Description:
 *      Read VCFs and output matrix file
 *  Arguments:
 *
 *  Returns:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-09  Jason Bacon Begin
 ***************************************************************************/

void    build_matrix(file_list_t *file_list, char *matrix_file)

{
    size_t  c;
    
    for (c = 0; c < file_list->count; ++c)
	puts(file_list->filename[c]);
    
    /*
     *  Read a call from each input file, output all those with the lowest
     *  chromosome/position or a . if the sample does not have a call there.
     *  Read the next call for each sample that was output and repeat until
     *  EOF on all files.
     */
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s filename-with-list-of-VCFs matrix-output-filename\n", argv[0]);
    exit(EX_USAGE);
}
