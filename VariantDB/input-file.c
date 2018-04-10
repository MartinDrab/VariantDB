
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "err.h"
#include "utils.h"
#include "file-utils.h"
#include "options.h"
#include "gen_dym_array.h"
#include "reads.h"
#include "input-file.h"


/************************************************************************/
/*                        HELPER FUNCTIONS                              */
/************************************************************************/

UTILS_TYPED_CALLOC_FUNCTION(ONE_READ)
UTILS_TYPED_CALLOC_FUNCTION(ACTIVE_REGION)


static const char *_read_line(const char *LineStart)
{
	while (*LineStart != '\n' && *LineStart != '\r' && *LineStart != 26 && *LineStart != '\0')
		++LineStart;

	return LineStart;
}


static const char *_advance_to_next_line(const char *LineEnd)
{
	while (*LineEnd == '\n' || *LineEnd == '\r')
		++LineEnd;

	return LineEnd;
}

typedef const char * cchar;


static boolean _fasta_read_seq_raw(char *Start, size_t Length, char **SeqStart, char **SeqEnd, cchar *Description, size_t *DescriptionLength)
{
	boolean ret = FALSE;

	ret = *Start == '>';
	if (ret) {
		size_t descrLen = 1;
		const char *descrStart = Start;

		while (Length > 0 && *Start != '\n' && *Start != '\r' && *Start != 26) {
			++Start;
			--Length;
			++descrLen;
		}

		*Description = descrStart;
		*DescriptionLength = descrLen - ((Length > 0) ? 1 : 0);
		if (Length > 0) {
			while (Length > 0 && (*Start == '\n' || *Start == '\r')) {
				++Start;
				--Length;
			}

			if (Length > 0 && *Start != 26) {
				*SeqStart = Start;
				while (Length > 0) {
					if (*Start == '>' && (*(Start - 1) == '\n' || *(Start - 1) == '\r'))
						break;

					++Start;
					--Length;
				}

				*SeqEnd = Start;
			} else {
				*SeqStart = NULL;
				*SeqEnd = NULL;
			}
		} else { 
			*SeqStart = NULL;
			*SeqEnd = NULL;
		}
	}

	return ret;
}



static ERR_VALUE _fasta_parse_description(const char *Description, size_t Length, char **Name, uint64_t *Pos)
{
	uint64_t tmpPos = 1;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const char *nameStart = NULL;
	size_t nameLen = 0;
	char *tmpName = NULL;

	*Name = NULL;
	*Pos = 1;
	ret = ERR_SUCCESS;
	assert(*Description == '>');
	++Description;
	--Length;
	while (Length > 0 && (*Description == ' ' || *Description == '\t')) {
		++Description;
		--Length;
	}

	if (Length > 0) {
		nameStart = Description;
		nameLen = 0;
		while (Length > 0 && *Description != ':') {
			++Description;
			--Length;
			++nameLen;
		}

		ret = utils_malloc((nameLen + 1)*sizeof(char), &tmpName);
		if (ret == ERR_SUCCESS) {
			memcpy(tmpName, nameStart, nameLen * sizeof(char));
			tmpName[nameLen] = '\0';
			if (Length > 1) {
				++Description;
				--Length;
				tmpPos = 0;
				while (Length > 0 && isdigit(*Description) != 0) {
					tmpPos = tmpPos * 10 + (*Description - '0');
					++Description;
					--Length;
				}
			}

			if (ret == ERR_SUCCESS) {
				*Pos = tmpPos;
				*Name = tmpName;
			}

			if (ret != ERR_SUCCESS)
				utils_free(tmpName);
		}
	}

	return ret;
}


static ERR_VALUE _fasta_read_seq(char *Start, size_t Length, char **NewStart, char **Seq, size_t *SeqLen, char **Name, uint64_t *Pos)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	char *seqStart = NULL;
	char *seqEnd = NULL;
	char *tmpSeq = NULL;
	size_t tmpSeqLen = 0;
	char *descr = NULL;
	size_t descrLen = 0;

	if (_fasta_read_seq_raw(Start, Length, &seqStart, &seqEnd, &descr, &descrLen)) {
		*Seq = NULL;
		*SeqLen = 0;
		if (seqStart != NULL) {
			ret = utils_malloc(seqEnd - seqStart + sizeof(char), &tmpSeq);
			if (ret == ERR_SUCCESS) {
				char *tmp = tmpSeq;

				while (seqStart != seqEnd) {
					if (*seqStart != '\n' && *seqStart != '\r' && *seqStart != 26) {
						*tmp = *seqStart;
						++tmp;
						++tmpSeqLen;
					}

					++seqStart;
				}

				tmpSeq[tmpSeqLen] = '\0';
				*NewStart = seqEnd;
				*Seq = tmpSeq;
				*SeqLen = tmpSeqLen;
				_fasta_parse_description(descr, descrLen, Name, Pos);
			}
		} else ret = ERR_NO_MORE_ENTRIES;
	} else ret = ERR_NO_MORE_ENTRIES;

	return ret;
}

/************************************************************************/
/*                        PUBLIC FUNCTIONS                              */
/************************************************************************/

ERR_VALUE fasta_load(const char *FileName, PFASTA_FILE FastaRecord)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_file_read(FileName, &FastaRecord->FileData, &FastaRecord->DataLength);
	if (ret == ERR_SUCCESS)
		FastaRecord->CurrentPointer = FastaRecord->FileData;

	return ret;
}



ERR_VALUE fasta_read_seq(PFASTA_FILE FastaRecord, PREFSEQ_DATA Data)
{
	size_t tmpLength = 0;
	char *tmpSeq = NULL;
	char *tmpName = NULL;
	uint64_t tmpPos = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	Data->StartPos = 0;
	Data->Name = NULL;
	ret = _fasta_read_seq(FastaRecord->CurrentPointer, FastaRecord->DataLength - (FastaRecord->CurrentPointer - FastaRecord->FileData), &FastaRecord->CurrentPointer, &tmpSeq, &tmpLength, &tmpName, &tmpPos);
	if (ret == ERR_SUCCESS) {
		Data->Sequence = tmpSeq;
		Data->Length = tmpLength;
		Data->Name = tmpName;
		Data->StartPos = tmpPos - 1;
	}

	return ret;
}


void fasta_free_seq(PREFSEQ_DATA Data)
{
	if (Data->Name != NULL)
		utils_free((char *)Data->Name);

	if (Data->Sequence != NULL)
		utils_free((char *)Data->Sequence);

	return;
}


void fasta_free(PFASTA_FILE FastaRecord)
{
	if (FastaRecord->DataLength > 0)
		utils_free(FastaRecord->FileData);

	return;
}


ERR_VALUE input_get_reads(const char *Filename, const CONFIDENT_REGION *Region, INPUT_READ_CALLBACK *Callback, void *Context)
{
	FILE *f = NULL;
	char line[4096];
	ONE_READ oneRead;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_fopen(Filename, FOPEN_MODE_READ, &f);
	if (ret == ERR_SUCCESS) {
		while (ret == ERR_SUCCESS && !feof(f) && !ferror(f)) {
			ret = utils_file_read_line(f, line, sizeof(line));
			if (ret == ERR_SUCCESS && *line != '@' && *line != '\0') {
				ret = read_create_from_sam_line(line, &oneRead);
				if (ret == ERR_SUCCESS &&
					((Region == NULL || strcmp(Region->Chrom, oneRead.Extension->RName) == 0) && Region->Start <= oneRead.Pos && oneRead.Pos < Region->End)) {
					if (!(oneRead.PosQuality < 20 || oneRead.Pos == (uint64_t)-1 ||
						oneRead.Extension->Flags.Bits.Unmapped ||
						oneRead.Extension->Flags.Bits.Supplementary ||
						oneRead.Extension->Flags.Bits.Duplicate ||
						oneRead.Extension->Flags.Bits.SecondaryAlignment)) {
						read_adjust(&oneRead, Region->Start, Region->End - Region->Start);
						ret = Callback(&oneRead, Context);
					}

					_read_destroy_structure(&oneRead);
				}
			}
		}

		utils_fclose(f);
	}
	
	return ret;
}


ERR_VALUE input_refseq_to_regions(const char *RefSeq, const size_t RefSeqLen, PACTIVE_REGION *Regions, size_t *Count)
{
	const char *regStart = NULL;
	const char *regEnd = NULL;
	size_t tmpArrayLen = 0;
	uint64_t numberOfNBases = 0;
	PACTIVE_REGION tmpArray = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	*Regions = NULL;
	*Count = 0;
	ret = ERR_SUCCESS;
	regStart = RefSeq;
	regEnd = RefSeq;
	for (size_t i = 0; i < RefSeqLen; ++i) {
		switch (*regEnd) {
			case 'A':
			case 'C':
			case 'G':
			case 'T':
				if (numberOfNBases > 0) {
					PACTIVE_REGION old = tmpArray;

					ret = utils_calloc_ACTIVE_REGION(tmpArrayLen + 1, &tmpArray);
					if (ret == ERR_SUCCESS) {
						PACTIVE_REGION newReg = tmpArray + tmpArrayLen;

						memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
						newReg->Sequence = regStart;
						newReg->Offset = (newReg->Sequence - RefSeq);
						newReg->Length = (regEnd - regStart);
						newReg->Type = artUnknown;
						++tmpArrayLen;
						if (old != NULL)
							utils_free(old);
					}

					numberOfNBases = 0;
					regStart = regEnd;
				}

				++regEnd;
				break;
			case 'N':
			case 'M':
			case 'R':
			case 'Y':
			case 'W':
			case 'S':
			case 'K':
			case 'V':
			case 'H':
			case 'D':
			case 'B':
				if (numberOfNBases == 0 && regEnd != regStart) {
					PACTIVE_REGION old = tmpArray;

					ret = utils_calloc_ACTIVE_REGION(tmpArrayLen + 1, &tmpArray);
					if (ret == ERR_SUCCESS) {
						PACTIVE_REGION newReg = tmpArray + tmpArrayLen;

						memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
						newReg->Sequence = regStart;
						newReg->Offset = (newReg->Sequence - RefSeq);
						newReg->Length = (regEnd - regStart);
						newReg->Type = artValid;
						++tmpArrayLen;
						if (old != NULL)
							utils_free(old);
					}
				}

				if (numberOfNBases == 0)
					regStart = regEnd;

				++regEnd;
				++numberOfNBases;
				break;
			default:
				printf("Unknown character in the reference sequence: %u\n", *regEnd);
				break;
		}

		if (ret != ERR_SUCCESS) {
			if (tmpArray != NULL)
				utils_free(tmpArray);

			break;
		}
	}

	if (ret == ERR_SUCCESS) {
		if (regEnd != regStart) {
			PACTIVE_REGION old = tmpArray;

			ret = utils_calloc_ACTIVE_REGION(tmpArrayLen + 1, &tmpArray);
			if (ret == ERR_SUCCESS) {
				PACTIVE_REGION newReg = tmpArray + tmpArrayLen;
				
				memcpy(tmpArray, old, tmpArrayLen*sizeof(ACTIVE_REGION));
				newReg->Sequence = regStart;
				newReg->Offset = (newReg->Sequence - RefSeq);
				newReg->Length = (regEnd - regStart);
				newReg->Type = (numberOfNBases == 0) ? artValid : artUnknown;
				tmpArrayLen++;
				if (old != NULL)
					utils_free(old);
			}
		}

		if (ret == ERR_SUCCESS) {
			*Regions = tmpArray;
			*Count = tmpArrayLen;
		}

		if (ret != ERR_SUCCESS) {
			if (tmpArray != NULL)
				utils_free(tmpArray);
		}
	}

	return ret;
}


ERR_VALUE input_get_region_by_offset(const PACTIVE_REGION Regions, const size_t Count, const uint64_t Offset, size_t *Index, uint64_t *RegionOffset)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	PACTIVE_REGION cur = Regions;
	uint64_t o = Offset;

	if (Offset >= cur->Offset) {
		ret = ERR_OFFSET_TOO_HIGH;
		for (size_t i = 0; i < Count; ++i) {
			if (cur->Length > o) {
				*RegionOffset = o;
				*Index = i;
				ret = ERR_SUCCESS;
				break;
			} else o -= cur->Length;

			++cur;
		}
	} else ret = ERR_OFFSET_TOO_LOW;

	return ret;
}


void input_free_regions(PACTIVE_REGION Regions, const size_t Count)
{
	if (Count > 0)
		utils_free(Regions);

	return;
}


static int _variant_comparator(const VCF_VARIANT *A, const VCF_VARIANT *B)
{
	int ret = 0;

	ret = strcmp(A->Chrom, B->Chrom);
	if (ret == 0)
		ret = (int)(A->Pos - B->Pos);

	return ret;
}


ERR_VALUE input_variant_create(const char *Chrom, const char *ID, unsigned long long Pos, const char *Ref, const char *Alt, unsigned long Quality, PVCF_VARIANT Variant)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	const size_t refLen = strlen(Ref);
	const size_t altLen = strlen(Alt);

	memset(Variant, 0, sizeof(VCF_VARIANT));
	ret = utils_copy_string(Chrom, &Variant->Chrom);
	if (ret == ERR_SUCCESS)
		ret = utils_copy_string(ID, &Variant->ID);

	if (ret == ERR_SUCCESS)
		ret = utils_copy_string(Ref, &Variant->Ref);

	if (ret == ERR_SUCCESS)
		ret = utils_copy_string(Alt, &Variant->Alt);

	if (ret == ERR_SUCCESS) {
		Variant->Pos = Pos;
		Variant->Quality = Quality;
		if (refLen == 1 && altLen == 1)
			Variant->Type = vcfvtSNP;
		else if (refLen > altLen && altLen == 1)
			Variant->Type = vcfvtDeletion;
		else if (altLen > refLen && refLen == 1)
			Variant->Type = vcfvtInsertion;
		else Variant->Type = vcfvtReplace;
	}

	if (ret != ERR_SUCCESS) {
		if (Variant->Alt != NULL)
			utils_free(Variant->Alt);

		if (Variant->Ref != NULL)
			utils_free(Variant->Ref);

		if (Variant->ID != NULL)
			utils_free(Variant->ID);

		if (Variant->Chrom != NULL)
			utils_free(Variant->Chrom);
	}

	return ret;
}

ERR_VALUE input_get_variants(const char *FileName, const VCF_VARIANT_FILTER *Filter, PGEN_ARRAY_VCF_VARIANT Array)
{
	size_t altLen = 0;
	char line[4096];
	FILE *f = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	POINTER_ARRAY_char fields;
	VCF_VARIANT v;
	unsigned long long pos = 0;
	unsigned long quality = 0;

	pointer_array_init_char(&fields, 140);
	ret = utils_fopen(FileName, FOPEN_MODE_READ, &f);
	if (ret == ERR_SUCCESS) {
		while (ret == ERR_SUCCESS && !feof(f) && !ferror(f)) {
			ret = utils_file_read_line(f, line, sizeof(line));
			if (ret == ERR_SUCCESS && *line != '\0' && *line != '#') {
				ret = utils_split(line, '\t', &fields);
				if (ret == ERR_SUCCESS) {
					pos = strtoull(fields.Data[1], NULL, 0) - 1;
					quality = strtoul(fields.Data[5], NULL, 0);
					
					char *alt = fields.Data[4];
					altLen = strlen(alt);
					ret = input_variant_create(fields.Data[0], fields.Data[2], pos, fields.Data[3], alt, quality, &v);
					if (Filter == NULL || input_variant_in_filter(Filter, &v))
						ret = dym_array_push_back_VCF_VARIANT(Array, v);

					if (ret == ERR_SUCCESS) {
						for (size_t i = 0; i < altLen; ++i) {
							if (alt[i] == ',') {
								alt[i] = '\0';
								ret = input_variant_create(fields.Data[0], fields.Data[2], pos, fields.Data[3], alt + i + 1, quality, &v);
								if (Filter == NULL || input_variant_in_filter(Filter, &v))
									ret = dym_array_push_back_VCF_VARIANT(Array, v);
							
								if (ret != ERR_SUCCESS)
									break;
							}
						}
					}

					utils_split_free(&fields);
				}
			}
		}

		utils_fclose(f);

	}

	pointer_array_finit_char(&fields);
	if (ret == ERR_SUCCESS)
		qsort(Array->Data, pointer_array_size(Array), sizeof(VCF_VARIANT), _variant_comparator);

	return ret;
}


void input_free_variant(const VCF_VARIANT *Variant)
{
	utils_free(Variant->Chrom);
	utils_free(Variant->ID);
	utils_free(Variant->Ref);
	utils_free(Variant->Alt);

	return;
}

void input_Free_variants(PGEN_ARRAY_VCF_VARIANT Array)
{
	const VCF_VARIANT *v = NULL;

	v = Array->Data;
	for (size_t i = 0; i < pointer_array_size(Array); ++i) {
		input_free_variant(v);
		++v;
	}

	dym_array_clear_VCF_VARIANT(Array);

	return;
}


boolean input_variant_in_filter(const VCF_VARIANT_FILTER *Filter, const VCF_VARIANT *Variant)
{
	long long int cmpResult = 0;
	boolean ret = FALSE;
	size_t intervalStart = 0;
	size_t intervalSize = Filter->RegionCount;
	size_t index = intervalStart + intervalSize / 2;

	ret = (intervalSize == 0);
	while (intervalSize > 0) {
		cmpResult = strcmp(Variant->Chrom, Filter->Regions[index].Chrom);
		if (cmpResult == 0) {
			cmpResult = (long long int)(Variant->Pos - Filter->Regions[index].Start);
			if (cmpResult >= 0) {
				ret = (Variant->Pos <= Filter->Regions[index].End);
				if (!ret) {
					intervalSize /= 2;
					if (cmpResult > 0)
						index += intervalSize;
					else index -= intervalSize;
				}
			} else {
				intervalSize /= 2;
				if (cmpResult > 0)
					index += intervalSize;
				else index -= intervalSize;
			}
		} else {
			intervalSize /= 2;
			if (cmpResult > 0)
				index += intervalSize;
			else index -= intervalSize;
		}
	}

	return ret;
}


boolean input_variant_normalize(const char *Reference, PVCF_VARIANT Variant)
{
	boolean ret = FALSE;

	switch (Variant->Type) {
		case vcfvtSNP:
			break;
		case vcfvtInsertion: {
			const size_t altLen = strlen(Variant->Alt);
			GEN_ARRAY_char altArray;
			const char *ref = Reference + Variant->Pos;

			dym_array_init_char(&altArray, 140);
			for (size_t i = 0; i < altLen; ++i)
				dym_array_push_back_char(&altArray, Variant->Alt[i]);

			while (*ref == altArray.Data[altArray.ValidLength - 1]) {
				--ref;
				memmove(altArray.Data + 1, altArray.Data, altArray.ValidLength - 1);
				altArray.Data[0] = *ref;
				Variant->Pos--;
			}

			memcpy(Variant->Alt, altArray.Data, altLen);
			Variant->Ref[0] = *ref;
			dym_array_finit_char(&altArray);
		} break;
		case vcfvtDeletion: {
			const size_t refLen = strlen(Variant->Ref);
			GEN_ARRAY_char refArray;
			const char *ref = Reference + Variant->Pos;

			dym_array_init_char(&refArray, 140);
			for (size_t i = 0; i < refLen; ++i)
				dym_array_push_back_char(&refArray, Variant->Ref[i]);

			while (*ref == refArray.Data[refArray.ValidLength - 1]) {
				--ref;
				memmove(refArray.Data + 1, refArray.Data, refArray.ValidLength - 1);
				refArray.Data[0] = *ref;
				Variant->Pos--;
			}

			memcpy(Variant->Ref, refArray.Data, refLen);
			Variant->Alt[0] = *ref;
			dym_array_finit_char(&refArray);
		} break;
		case vcfvtReplace:
			break;
		default:
			break;
	}

	return ret;
}


static int _bed_comparator(const CONFIDENT_REGION *A, const CONFIDENT_REGION *B)
{
	int ret = 0;

	ret = strcmp(A->Chrom, B->Chrom);
	if (ret == 0)
		ret = (int)(A->Start - B->Start);

	return ret;
}


ERR_VALUE input_get_bed(const char *FileName, const char *Chrom, PGEN_ARRAY_CONFIDENT_REGION Array)
{
	char line[4096];
	FILE *f = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;
	POINTER_ARRAY_char fields;
	CONFIDENT_REGION cr;

	pointer_array_init_char(&fields, 140);
	ret = utils_fopen(FileName, FOPEN_MODE_READ, &f);
	if (ret == ERR_SUCCESS) {
		while (ret == ERR_SUCCESS && !feof(f) && !ferror(f)) {
			ret = utils_file_read_line(f, line, sizeof(line));
			if (ret == ERR_SUCCESS && *line != '\0' && *line != '#') {
				ret = utils_split(line, '\t', &fields);
				if (ret == ERR_SUCCESS) {
					cr.Chrom = fields.Data[0];
					cr.Start = strtoull(fields.Data[1], NULL, 0);
					cr.End = strtoull(fields.Data[2], NULL, 0);
					if (Chrom == NULL || *Chrom == '\0' || strcmp(cr.Chrom, Chrom) == 0) {
						ret = dym_array_push_back_CONFIDENT_REGION(Array, cr);
						if (ret == ERR_SUCCESS)
							fields.Data[0] = NULL;
					}

					utils_split_free(&fields);
				}
			}
		}

		utils_fclose(f);

	}

	pointer_array_finit_char(&fields);
	if (ret == ERR_SUCCESS)
		qsort(Array->Data, pointer_array_size(Array), sizeof(CONFIDENT_REGION), _bed_comparator);
	
	return ret;
}


void input_free_bed(PGEN_ARRAY_CONFIDENT_REGION Array)
{
	const CONFIDENT_REGION *cr = NULL;

	cr = Array->Data;
	for (size_t i = 0; i < pointer_array_size(Array); ++i) {
		utils_free(cr->Chrom);
		++cr;
	}

	dym_array_clear_CONFIDENT_REGION(Array);

	return;
}
