

#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "input-file.h"
#include "variantdb.h"


static char *_refFile = NULL;
static char *_samFile = NULL;
static char *_vcfFile = NULL;
static char *_bedFile = NULL;
static char *_chromosome = NULL;
static uint64_t _regionStart = 0;
static uint64_t _regionEnd = (uint64_t)-1;
static boolean _help = FALSE;
static boolean _verbose = FALSE;


static void _cmd_option_init(void)
{
	CMD_OPTION_INIT(VDB_OPTION_REF_FILE, String, "");
	CMD_OPTION_INIT(VDB_OPTION_SAM_FILE, String, "");
	CMD_OPTION_INIT(VDB_OPTION_VCF_FILE, String, "");
	CMD_OPTION_INIT(VDB_OPTION_BED_FILE, String, "");
	CMD_OPTION_INIT(VDB_OPTION_CHROM, String, "1");
	CMD_OPTION_INIT(VDB_OPTION_START, UInt64, 0);
	CMD_OPTION_INIT(VDB_OPTION_STOP, UInt64, (uint64_t)-1);
	CMD_OPTION_INIT(VDB_OPTION_HELP, Boolean, FALSE);
	CMD_OPTION_INIT(VDB_OPTION_VERBOSE, Boolean, FALSE);

	return;
}


static ERR_VALUE _cmd_optiion_parse(void)
{
	CMD_OPTION_GET(VDB_OPTION_REF_FILE, String, &_refFile);
	CMD_OPTION_GET(VDB_OPTION_SAM_FILE, String, &_samFile);
	CMD_OPTION_GET(VDB_OPTION_VCF_FILE, String, &_vcfFile);
	CMD_OPTION_GET(VDB_OPTION_BED_FILE, String, &_bedFile);
	CMD_OPTION_GET(VDB_OPTION_CHROM, String, &_chromosome);
	CMD_OPTION_GET(VDB_OPTION_START, UInt64, &_regionStart);
	CMD_OPTION_GET(VDB_OPTION_STOP, UInt64, &_regionEnd);
	CMD_OPTION_GET(VDB_OPTION_HELP, Boolean, &_help);
	CMD_OPTION_GET(VDB_OPTION_VERBOSE, Boolean, &_verbose);

	if (_help)
		return ERR_SUCCESS;

	if (*_refFile == '\0') {
		fprintf(stderr, "[ERROR]: The reference sequence file was not specified (--%s)\n", VDB_OPTION_REF_FILE);
		return ERR_INTERNAL_ERROR;
	}

	if (*_samFile == '\0') {
		fprintf(stderr, "[ERROR]: The SAM file was not specified (--%s)\n", VDB_OPTION_SAM_FILE);
		return ERR_INTERNAL_ERROR;
	}

	if (*_vcfFile == '\0') {
		fprintf(stderr, "[ERROR]: The VCF file was not specified (--%s)\n", VDB_OPTION_VCF_FILE);
		return ERR_INTERNAL_ERROR;
	}

	if (*_bedFile == '\0')
		fprintf(stderr, "[WARNING]: The BED file was not specified (--%s). Treating all regions as confident\n", VDB_OPTION_BED_FILE);

	if (_regionStart >= _regionEnd) {
		fprintf(stderr, "[ERROR]: The specified region (--%s, --%s) is not an interval\n", VDB_OPTION_START, VDB_OPTION_STOP);
		return ERR_INTERNAL_ERROR;
	}

	return ERR_SUCCESS;
}


static FASTA_FILE refFile;
static REFSEQ_DATA refData;
static boolean _fastaLoaded = FALSE;



int main(int argc, char **argv)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_allocator_init(1);
	if (ret == ERR_SUCCESS) {
		ret = options_module_init(37);
		if (ret == ERR_SUCCESS) {
			_cmd_option_init();
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS)
				ret = _cmd_optiion_parse();
			
			if (ret == ERR_SUCCESS) {
				ret = fasta_load(_refFile, &refFile);
				if (ret == ERR_SUCCESS) {
					ret = fasta_read_seq(&refFile, &refData);
					if (ret != ERR_SUCCESS)
						fasta_free(&refFile);
				}
			}
			
			if (ret == ERR_SUCCESS) {
				_fastaLoaded = TRUE;
			}

			if (_fastaLoaded) {
				fasta_free_seq(&refData);
				fasta_free(&refFile);
			}

			options_module_finit();
		}
	}

	return 0;
}
