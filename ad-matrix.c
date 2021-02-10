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
#include <stdbool.h>

#include <tsvio.h>
#include <vcfio.h>
#include <biostring.h>
#include <errno.h>

#include "ad-matrix.h"

int     main(int argc,char *argv[])

{
    file_list_t file_list;
    char        *list_filename = argv[1],
		*matrix_filename_stem = argv[2];
    
    if ( argc != 3 )
	usage(argv);
    open_files(list_filename, &file_list, "r");
    build_matrix(&file_list, matrix_filename_stem);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Read a list of VCF files from filename and open all files with
 *      the given fopen() mode
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-09  Jason Bacon Begin
 ***************************************************************************/

void    open_files(char *list_filename, file_list_t *file_list, char *mode)

{
    FILE        *fp;
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
    puts("All files opened.");
}


/***************************************************************************
 *  Description:
 *      Read VCFs and output matrix file
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-09  Jason Bacon Begin
 ***************************************************************************/

void    build_matrix(file_list_t *file_list, char *matrix_stem)

{
    size_t      c,
		low_pos,
		open_count,
		rows = 0;
    vcf_call_t  *vcf_call;
    int         chr_cmp;
    char        *ref_count,
		*alt_count,
		*ref_alt_count,
		*low_chrom,
		ref_matrix_pipe[PATH_MAX + 1],
		ref_alt_matrix_pipe[PATH_MAX + 1];
    FILE        *ref_matrix_fp,
		*ref_alt_matrix_fp;
    
    vcf_call = (vcf_call_t *)malloc(file_list->count * sizeof(vcf_call_t));
    if ( vcf_call == NULL )
    {
	fprintf(stderr, "build_matrix(): Could not allocate vcf_call array.\n");
	fprintf(stderr, "Size = %zu\n", file_list->count * sizeof(vcf_call_t));
	exit(EX_UNAVAILABLE);
    }
    
    /* Use a lower compression level than default 6 so xz can keep up */
    snprintf(ref_matrix_pipe, PATH_MAX, "xz -4 - > %s-ref.tsv.xz", matrix_stem);
    if ( (ref_matrix_fp = popen(ref_matrix_pipe, "w")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", ref_matrix_pipe,
		strerror(errno));
	exit(EX_CANTCREAT);
    }
    
    snprintf(ref_alt_matrix_pipe, PATH_MAX, "xz -4 - > %s-ref+alt.tsv.xz", matrix_stem);
    if ( (ref_alt_matrix_fp = popen(ref_alt_matrix_pipe, "w")) == NULL )
    {
	fprintf(stderr, "Cannot open %s: %s\n", ref_alt_matrix_pipe,
		strerror(errno));
	exit(EX_CANTCREAT);
    }
    
    /*
     *  Read a call from each input file, output all those with the lowest
     *  chromosome/position or a . if the sample does not have a call there.
     *  Read the next call for each sample that was output and repeat until
     *  EOF on all files.
     */

    /* First call from each sample file */
    puts("Reading first call from each sample...");
    for (c = 0; c < file_list->count; ++c)
    {
	vcf_call_init(&vcf_call[c], 16, 32, 64);
	if ( vcf_read_ss_call(file_list->fp[c], &vcf_call[c]) == VCF_OK )
	{
#ifdef DEBUG
	    fprintf(ref_matrix_fp, "%zu %s %s %zu %s\n",
		    c, file_list->filename[c],
		    VCF_CHROMOSOME(&vcf_call[c]),
		    VCF_POS(&vcf_call[c]), VCF_SINGLE_SAMPLE(&vcf_call[c]));
	    fprintf(ref_alt_matrix_fp, "%zu %s %s %zu %s\n",
		    c, file_list->filename[c],
		    VCF_CHROMOSOME(&vcf_call[c]),
		    VCF_POS(&vcf_call[c]), VCF_SINGLE_SAMPLE(&vcf_call[c]));
#endif
	}
	else
	{
	    fprintf(stderr,
		    "build_matrix(): Failed to read VCF call from %s.\n",
		    file_list->filename[c]);
	    exit(EX_DATAERR);
	}
    }
    puts("First calls read.");

    open_count = file_list->count;
    while ( open_count > 0 )
    {
	/*
	 *  Find lowest pos among all samples
	 */
	
	/* Skip over finished sample files */
	for (c = 0; file_list->fp[c] == NULL; ++c)
	    ;
	
	/* Assume first sample has lowest position than scan the rest */
	low_pos = VCF_POS(&vcf_call[c]);
	low_chrom = VCF_CHROMOSOME(&vcf_call[c]);
	for (c = c + 1; c < file_list->count; ++c)
	{
	    chr_cmp = chromosome_name_cmp(VCF_CHROMOSOME(&vcf_call[c]),low_chrom);
	    if ( (file_list->fp[c] != NULL) && ((chr_cmp < 0) ||
		    ((chr_cmp == 0) && (VCF_POS(&vcf_call[c]) < low_pos))) )
	    {
		low_pos = VCF_POS(&vcf_call[c]);
		low_chrom = VCF_CHROMOSOME(&vcf_call[c]);
	    }
	}
	
	/* Output row for low pos, read next call for represented samples */
	fprintf(ref_matrix_fp, "%s\t%zu\t", low_chrom, low_pos);
	fprintf(ref_alt_matrix_fp, "%s\t%zu\t", low_chrom, low_pos);
	for (c = 0; c < file_list->count; ++c)
	{
	    if ( (file_list->fp[c] != NULL) &&
		 (VCF_POS(&vcf_call[c]) == low_pos) )
	    {
		/* Find start of ref count */
		ref_count = VCF_SINGLE_SAMPLE(&vcf_call[c]);
		strsep(&ref_count, ":");
		
		/* null-terminate ref count */
		alt_count = ref_count;
		strsep(&alt_count, ",");
		
		/* Find start of depth, already null-terminated, last field */
		ref_alt_count = alt_count;
		strsep(&ref_alt_count, ":");
		fprintf(ref_matrix_fp, "%s\t", ref_count);
		fprintf(ref_alt_matrix_fp, "%s\t", ref_alt_count);
		if ( vcf_read_ss_call(file_list->fp[c], &vcf_call[c]) ==
			VCF_READ_EOF )
		{
		    fprintf(stderr, "Closing %zu %s\n", c,
			    file_list->filename[c]);
		    fclose(file_list->fp[c]);
		    file_list->fp[c] = NULL;
		    --open_count;
		}
	    }
	    else
	    {
		fprintf(ref_matrix_fp, ".\t");
		fprintf(ref_alt_matrix_fp, ".\t");
	    }
	}
	putc('\n', ref_matrix_fp);
	putc('\n', ref_alt_matrix_fp);
	
#ifdef DEBUG
	for (c = 0; c < file_list->count; ++c)
	{
	    if ( file_list->fp[c] != NULL )
	    {
		fprintf(ref_matrix_fp, "%zu %s %s %zu %s\n",
			c, file_list->filename[c],
			VCF_CHROMOSOME(&vcf_call[c]),
			VCF_POS(&vcf_call[c]),
			VCF_SINGLE_SAMPLE(&vcf_call[c]));
		fprintf(ref_alt_matrix_fp, "%zu %s %s %zu %s\n",
			c, file_list->filename[c],
			VCF_CHROMOSOME(&vcf_call[c]),
			VCF_POS(&vcf_call[c]),
			VCF_SINGLE_SAMPLE(&vcf_call[c]));
	    }
	    else
	    {
		fprintf(ref_matrix_fp, "%zu %s EOF\n",
			c, file_list->filename[c]);
		fprintf(ref_alt_matrix_fp, "%zu %s EOF\n",
			c, file_list->filename[c]);
	    }
	}
#endif

	if ( ++rows % 1000 == 0 )
	    fprintf(stderr, "%zu\r", rows);
    }
    pclose(ref_matrix_fp);
    pclose(ref_alt_matrix_fp);
    fprintf(stderr, "Done!\n");
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s filename-with-list-of-VCFs matrix-output-stem\n", argv[0]);
    fprintf(stderr, "Two matrix files are produced, named\n");
    fprintf(stderr, "<matrix-output-stem>-ref.tsv and <matrix-output-stem>-ref+alt.tsv\n");
    exit(EX_USAGE);
}
