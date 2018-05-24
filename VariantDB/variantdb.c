

#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "options.h"
#include "khash.h"
#include "input-file.h"
#include "ssw.h"
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
static VCF_VARIANT_FILTER variantFilter;
static CONFIDENT_REGION region;
static GEN_ARRAY_VCF_VARIANT variants;
static boolean _variantsLoaded = FALSE;
static GEN_ARRAY_CONFIDENT_REGION confidentRegions;
static boolean _bedLoaded = FALSE;


KHASH_MAP_INIT_INT64(VariantTableType, PVCF_VARIANT);

khash_t(VariantTableType) *_variantTable = NULL;
khash_t(VariantTableType) *_nearVariantTable = NULL;

static size_t _readsProcessed = 0;

static ERR_VALUE _on_read_callback(const ONE_READ *Read, void *Context)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const char *ref = refData.Sequence + Read->Pos - refData.StartPos;
	char *opString = NULL;
	size_t opStringSize = 0;
	unsigned long long currentPos = Read->Pos;
	unsigned long long variantPos = 0;
	const char *currentOp = NULL;
	GEN_ARRAY_char refArray;
	GEN_ARRAY_char altArray;
	size_t readSeqIndex = 0;
	PVCF_VARIANT v = NULL;
	khiter_t it;

	dym_array_init_char(&refArray, 140);
	dym_array_init_char(&altArray, 140);
	ret = ssw_clever(ref, Read->ReadSequenceLen, Read->ReadSequence, Read->ReadSequenceLen, 2, -1, -1, &opString, &opStringSize);
	if (ret == ERR_SUCCESS) {
		currentOp = opString;
		while (*currentOp != '\0') {
			switch (*currentOp) {
				case 'I':
					if (variantPos == 0)
						variantPos = currentPos;

					dym_array_push_back_char(&altArray, Read->ReadSequence[readSeqIndex]);
					++readSeqIndex;
					break;
				case 'D':
					if (variantPos == 0)
						variantPos = currentPos;

					dym_array_push_back_char(&refArray, *ref);
					++ref;
					++currentPos;
					it = kh_get(VariantTableType, _variantTable, currentPos);
					if (it != kh_end(_variantTable))
						kh_value(_variantTable, it)->TotalReadsAtPosition++;
					break;
				case 'X':
					if (variantPos == 0)
						variantPos = currentPos;

					dym_array_push_back_char(&refArray, *ref);
					++ref;
					++currentPos;
					dym_array_push_back_char(&altArray, Read->ReadSequence[readSeqIndex]);
					++readSeqIndex;
					it = kh_get(VariantTableType, _variantTable, currentPos);
					if (it != kh_end(_variantTable))
						kh_value(_variantTable, it)->TotalReadsAtPosition++;
					break;
				case 'M':
					if (variantPos != 0) {
						if (gen_array_size(&altArray) == 0) {
							variantPos--;
							dym_array_push_back_char(&altArray, refData.Sequence[variantPos]);
							dym_array_push_back_char(&refArray, '\0');
							memmove(refArray.Data + 1, refArray.Data, refArray.ValidLength - 1);
							refArray.Data[0] = refData.Sequence[variantPos];
						}

						if (gen_array_size(&refArray) == 0) {
							variantPos--;
							dym_array_push_back_char(&refArray, refData.Sequence[variantPos]);
							dym_array_push_back_char(&altArray, '\0');
							memmove(altArray.Data + 1, altArray.Data, altArray.ValidLength - 1);
							altArray.Data[0] = refData.Sequence[variantPos];
						}

						dym_array_push_back_char(&refArray, '\0');
						dym_array_push_back_char(&altArray, '\0');
						ret = utils_malloc(sizeof(VCF_VARIANT), &v);
						if (ret == ERR_SUCCESS) {
							memset(v, 0, sizeof(VCF_VARIANT));
							ret = input_variant_create(Read->Extension->RName, NULL, variantPos, refArray.Data, altArray.Data, 30, v);
							if (ret == ERR_SUCCESS) {
								input_variant_normalize(refData.Sequence, v);
								it = kh_get(VariantTableType, _variantTable, v->Pos);
								if (it != kh_end(_variantTable) &&
									input_variant_equal(v, kh_value(_variantTable, it))) {
									kh_value(_variantTable, it)->ReadSupport++;
									input_free_variant(v);
									utils_free(v);
									v = NULL;
								} else {
									int res = 0;
									PVCF_VARIANT tmp = NULL;
									PVCF_VARIANT tmp2 = NULL;

									it = kh_put(VariantTableType, _nearVariantTable, v->Pos, &res);
									if (res == 0) {
										tmp = kh_value(_nearVariantTable, it);
										tmp2 = tmp;
										while (tmp2 != NULL) {
											if (input_variant_equal(tmp2, v)) {
												++tmp2->ReadSupport;
												input_free_variant(v);
												utils_free(v);
												v = NULL;
												break;
											}

											tmp2 = tmp2->Alternative;
										}

										if (tmp2 == NULL)
											v->Alternative = tmp;
									}

									if (v != NULL)
										kh_value(_nearVariantTable, it) = v;

								}
							}
						}

						dym_array_clear_char(&refArray);
						dym_array_clear_char(&altArray);
						variantPos = 0;
					}

					++readSeqIndex;
					++ref;
					++currentPos;
					it = kh_get(VariantTableType, _variantTable, currentPos);
					if (it != kh_end(_variantTable))
						kh_value(_variantTable, it)->TotalReadsAtPosition++;
					break;
			}

			++currentOp;
		}

		utils_free(opString);
	}

	dym_array_finit_char(&altArray);
	dym_array_finit_char(&refArray);
	_readsProcessed++;
	if (_readsProcessed % 10000 == 0) {
		fputc('.', stderr);
		fflush(stderr);
	}

	return ret;
}


int main(int argc, char **argv)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	_variantTable = kh_init(VariantTableType);
	_nearVariantTable = kh_init(VariantTableType);
	ret = utils_allocator_init(1);
	if (ret == ERR_SUCCESS) {
		ret = options_module_init(37);
		if (ret == ERR_SUCCESS) {
			_cmd_option_init();
			ret = options_parse_command_line(argc - 1, argv + 1);
			if (ret == ERR_SUCCESS)
				ret = _cmd_optiion_parse();

			if (ret == ERR_SUCCESS && !_help) {
				fprintf(stderr, "[INFO]: Loading the reference...\n");
				ret = fasta_load(_refFile, &refFile);
				if (ret == ERR_SUCCESS) {
					ret = fasta_read_seq(&refFile, &refData);
					if (ret != ERR_SUCCESS)
						fasta_free(&refFile);
				}

				fprintf(stderr, "[INFO]: Chromosome is set to %s\n", _chromosome);
				fprintf(stderr, "[INFO]: Area start is set to %llu\n", _regionStart);
				fprintf(stderr, "[INFO]: Area end is set to %llu\n", _regionEnd);
				region.Chrom = _chromosome;
				region.Start = _regionStart;
				region.End = _regionEnd;
				_fastaLoaded = (ret == ERR_SUCCESS);
				if (ret == ERR_SUCCESS && *_bedFile != '\0') {
					fprintf(stderr, "[INFO]: Loading the BED...\n");
					dym_array_init_CONFIDENT_REGION(&confidentRegions, 150);
					ret = input_get_bed(_bedFile, &region, &confidentRegions);
					_bedLoaded = (ret == ERR_SUCCESS);
					if (_bedLoaded)
						fprintf(stderr, "[INFO]: %zu confidence regions loaded\n", gen_array_size(&confidentRegions));
				}

				if (ret == ERR_SUCCESS) {
					fprintf(stderr, "[INFO]: Loading the VCF...\n");
					if (!_bedLoaded) {
						variantFilter.RegionCount = 1;
						variantFilter.Regions = &region;
					} else {
						variantFilter.RegionCount = gen_array_size(&confidentRegions);
						variantFilter.Regions = confidentRegions.Data;
					}

					dym_array_init_VCF_VARIANT(&variants, 150);
					ret = input_get_variants(_vcfFile, &variantFilter, &variants);
					if (ret == ERR_SUCCESS) {
						khiter_t it;
						int res = 0;
						PVCF_VARIANT v = variants.Data;
						PVCF_VARIANT tmp = NULL;

						for (size_t i = 0; i < gen_array_size(&variants); ++i) {
							if (input_variant_normalize(refData.Sequence, v))
								fprintf(stderr, "[ERROR]: Normalized:\t%llu\t%s\t%s\n", v->Pos + 1, v->Ref, v->Alt);
							
							it = kh_put(VariantTableType, _variantTable, v->Pos, &res);
							if (res == 0) {
								tmp = kh_value(_variantTable, it);
								v->Alternative = tmp;
							}

							kh_value(_variantTable, it) = v;
							++v;
						}
					}

					_variantsLoaded = (ret == ERR_SUCCESS);
					if (_variantsLoaded)
						fprintf(stderr, "[INFO]: %zu variants loaded\n", gen_array_size(&variants));
				}

				if (ret == ERR_SUCCESS) {
					fprintf(stderr, "[INFO]: Processing the reads...\n");
					ret = input_get_reads(_samFile, &region, _on_read_callback, NULL);
				}

				if (ret == ERR_SUCCESS) {
					PVCF_VARIANT v = variants.Data;
					PVCF_VARIANT tmp = NULL;

					fprintf(stderr, "\n");
					fprintf(stderr, "[INFO]: Processing variants...\n");
					for (size_t i = 0; i < gen_array_size(&variants); ++i) {
						fprintf(stdout, "%s\t%llu\t%s\t%s\t%s\t%zu\t%zu\n", v->Chrom, v->Pos + 1, v->ID, v->Ref, v->Alt, v->ReadSupport, v->TotalReadsAtPosition);
						for (size_t j = v->Pos - 10; j < v->Pos + 10; ++j) {
							khiter_t it = kh_get(VariantTableType, _nearVariantTable, j);
							if (it != kh_end(_nearVariantTable)) {
								tmp = kh_val(_nearVariantTable, it);
								fprintf(stdout, "\t%s\t%llu\t%s\t%s\t%s\t%zu\t%zu\n", tmp->Chrom, tmp->Pos + 1, tmp->ID, tmp->Ref, tmp->Alt, tmp->ReadSupport, tmp->TotalReadsAtPosition);
							}
						}
						++v;
					}
				}

				if (_bedLoaded) {
					fprintf(stderr, "[INFO]: Freeing the BED...\n");
					input_free_bed(&confidentRegions);
					dym_array_finit_CONFIDENT_REGION(&confidentRegions);
				}

				if (_variantsLoaded) {
					fprintf(stderr, "[INFO]: Freeing the VCF...\n");
					input_Free_variants(&variants);
					dym_array_finit_VCF_VARIANT(&variants);
				}

				if (_fastaLoaded) {
					fprintf(stderr, "[INFO]: Freeing the reference...\n");
					fasta_free_seq(&refData);
					fasta_free(&refFile);
				}
			} else if (_help)
				options_print_help();
			else fprintf(stderr, "[INFO]: Use variantdb -h for help\n");

			options_module_finit();
		}
	}

	if (ret != ERR_SUCCESS)
		fprintf(stderr, "[ERROR]: The operation failed with an error code %u\n", ret);

	return ret;
}
