
#include <stdlib.h>
#include <inttypes.h>
#include "err.h"
#include "utils.h"
#include "file-utils.h"
#include "reads.h"



/************************************************************************/
/*                          HELPER FUNCTIONS                            */
/************************************************************************/

UTILS_TYPED_CALLOC_FUNCTION(ONE_READ)


static const char *_sam_read_field(const char *Start)
{
	while (*Start != '\0' && *Start != '\r' && *Start != '\n' && *Start != '\t' && *Start != 26)
		++Start;

	return Start;
}


static const char *_sam_read_string_field(const char *Start, char **String, size_t *Length)
{
	char *tmpString = NULL;
	const char *end = NULL;
	size_t len = 0;

	end = _sam_read_field(Start);
	len = (end - Start);
	if (len > 0) {
		if (utils_malloc((len + 1)*sizeof(char), &tmpString) == ERR_SUCCESS) {
			memcpy(tmpString, Start, len*sizeof(char));
			tmpString[len] = '\0';
			if (String != NULL) {
				*String = tmpString;
				*Length = len;
			} else utils_free(tmpString);
		} else end = NULL;
	} else end = NULL;

	return end;
}


static const char *_sam_read_uint_field(const char *Start, uint32_t *Value)
{
	uint32_t tmpValue = 0;
	const char *end = NULL;
	size_t len = 0;

	end = _sam_read_field(Start);
	len = (end - Start);
	if (len > 0) {
		char *tmpEnd = NULL;

		tmpValue = (uint32_t)strtoul(Start, &tmpEnd, 0);
		if (end == tmpEnd) {
			if (Value != NULL)
				*Value = tmpValue;
		} else end = NULL;
	} else end = NULL;

	return end;
}


static const char *_sam_read_int_field(const char *Start, int32_t *Value)
{
	int32_t tmpValue = 0;
	const char *end = NULL;
	size_t len = 0;

	end = _sam_read_field(Start);
	len = (end - Start);
	if (len > 0) {
		char *tmpEnd = NULL;

		tmpValue = (int32_t)strtol(Start, &tmpEnd, 0);
		if (end == tmpEnd) {
			if (Value != NULL)
				*Value = tmpValue;
		} else end = NULL;
	}
	else end = NULL;

	return end;
}



void _read_destroy_structure(PONE_READ Read)
{
	Read->Quality -= Read->Offset;
	Read->ReadSequence -= Read->Offset;

	if (Read->Quality != NULL)
		utils_free(Read->Quality);

	if (Read->ReadSequence != NULL)
		utils_free(Read->ReadSequence);

	if (Read->Extension->RNext != NULL)
		utils_free(Read->Extension->RNext);

	if (Read->Extension->CIGAR != NULL)
		utils_free(Read->Extension->CIGAR);

	if (Read->Extension->RName != NULL)
		utils_free(Read->Extension->RName);

	if (Read->Extension->TemplateName != NULL)
		utils_free(Read->Extension->TemplateName);

	utils_free(Read->Extension);

	return;
}



/************************************************************************/
/*                          PUBLIC FUNCTIONS                            */
/************************************************************************/


void read_quality_encode(PONE_READ Read)
{
	for (size_t i = 0; i < Read->ReadSequenceLen; ++i)
		Read->Quality[i] += 33;

	return;
}


void read_quality_decode(PONE_READ Read)
{
	for (size_t i = 0; i < Read->ReadSequenceLen; ++i)
		Read->Quality[i] -= 33;

	return;
}


void read_write_fastq(FILE *Stream, const ONE_READ *Read)
{
	fprintf(Stream, "@%s\n", Read->Extension->TemplateName);
	fprintf(Stream, "%.*s\n", Read->ReadSequenceLen, Read->ReadSequence);
	fprintf(Stream, "+%s\n", Read->Extension->TemplateName);
	fprintf(Stream, "%.*s\n", Read->ReadSequenceLen, Read->Quality);

	return;
}


static void _advance_to_line_end(const char **pCursor)
{
	const char *tmp = *pCursor;

	while (*tmp != '\0' && *tmp != '\r' && *tmp != '\n' && *tmp != 26)
		++tmp;

	*pCursor = tmp;

	return;
}


static boolean _advance_to_next_line(const char **pCursor)
{
	boolean ret = TRUE;
	const char *tmp = *pCursor;

	while (*tmp == '\n' || *tmp == '\r')
		++tmp;

	*pCursor = tmp;
	ret = (*tmp != 26 && *tmp != '\0');

	return ret;
}


ERR_VALUE read_create_from_fastq(const char *Block, const char **NewBlock, PONE_READ Read)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Read, 0, sizeof(ONE_READ));
	ret = utils_malloc(sizeof(ONE_READ_EXTENSION), &Read->Extension);
	if (ret == ERR_SUCCESS) {
		memset(Read->Extension, 0, sizeof(ONE_READ_EXTENSION));
		if (*Block == '@') {
			size_t templateSize = 0;
			const char *lineEnd = Block + 1;

			_advance_to_line_end(&lineEnd);
			templateSize = lineEnd - Block - 1;
			ret = utils_calloc_char(templateSize + 1, &Read->Extension->TemplateName);
			if (ret == ERR_SUCCESS) {
				memcpy(Read->Extension->TemplateName, Block + 1, templateSize * sizeof(char));
				Read->Extension->TemplateName[templateSize] = '\0';
				if (_advance_to_next_line(&lineEnd)) {
					size_t seqLen = 0;

					Block = lineEnd;
					_advance_to_line_end(&lineEnd);
					seqLen = lineEnd - Block;
					ret = utils_calloc_char(seqLen + 1, &Read->ReadSequence);
					if (ret == ERR_SUCCESS) {
						memcpy(Read->ReadSequence, Block, seqLen * sizeof(char));
						Read->ReadSequence[seqLen] = '\0';
						Read->ReadSequenceLen = (uint32_t)seqLen;
						if (_advance_to_next_line(&lineEnd) && *lineEnd == '+') {
							Block = lineEnd;
							_advance_to_line_end(&lineEnd);
							if (_advance_to_next_line(&lineEnd)) {
								size_t qLen = 0;

								Block = lineEnd;
								_advance_to_line_end(&lineEnd);
								qLen = lineEnd - Block;
								if (qLen == Read->ReadSequenceLen) {
									ret = utils_calloc_char(Read->ReadSequenceLen + 1, &Read->Quality);
									if (ret == ERR_SUCCESS) {
										memcpy(Read->Quality, Block, Read->ReadSequenceLen * sizeof(char));
										Read->Quality[Read->ReadSequenceLen] = '\0';
										read_quality_decode(Read);
										_advance_to_next_line(&lineEnd);
										*NewBlock = lineEnd;
									}
								} else ret = ERR_FASTQ_LEN_MISMATCH;
							} else ret = ERR_FASTQ_NO_QUALITY;
						} else ret = ERR_FASTQ_NO_PLUS;

						if (ret != ERR_SUCCESS)
							utils_free(Read->ReadSequence);
					}
				} else ret = ERR_FASTQ_NO_SEQ;

				if (ret != ERR_SUCCESS)
					utils_free(Read->Extension->TemplateName);
			}
		} else ret = ERR_FASTQ_NO_DESCRIPTION;
		
		if (ret != ERR_SUCCESS)
			utils_free(Read->Extension);
	}

	return ret;
}

void read_write_sam(FILE *Stream, const ONE_READ *Read)
{
	fprintf(Stream, "%s\t%u\t%s\t%" PRId64 "\t%u\t%s\t%s\t%" PRId64 "\t%i\t%.*s\t%.*s\n", Read->Extension->TemplateName, Read->Extension->Flags.Value, Read->Extension->RName, Read->Pos + 1, Read->PosQuality, Read->Extension->CIGAR, Read->Extension->RNext, Read->Extension->PNext, Read->Extension->TLen, (int)Read->ReadSequenceLen, Read->ReadSequence, (int)Read->ReadSequenceLen, Read->Quality);

	return;
}


ERR_VALUE read_create_from_sam_line(const char *Line, PONE_READ Read)
{
	uint32_t tmp32;
	size_t tmpStringLen = 0;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	memset(Read, 0, sizeof(ONE_READ));
	ret = utils_malloc(sizeof(ONE_READ_EXTENSION), &Read->Extension);
	if (ret == ERR_SUCCESS) {
		Line = _sam_read_string_field(Line, &Read->Extension->TemplateName, &tmpStringLen);
		if (Line != NULL && *Line == '\t')
			++Line;
		else ret = ERR_SAM_INVALID_RNAME;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_uint_field(Line, &tmp32);
		if (Line != NULL && *Line == '\t') {
			++Line;
			Read->Extension->Flags.Value = tmp32;
		} else ret = ERR_SAM_INVALID_FLAG;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_string_field(Line, &Read->Extension->RName, &tmpStringLen);
		if (Line != NULL && *Line == '\t')
			++Line;
		else ret = ERR_SAM_INVALID_RNAME;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_uint_field(Line, &tmp32);
		if (Line != NULL && *Line == '\t') {
			++Line;
			Read->Pos = tmp32;
			Read->Pos--;
		} else ret = ERR_SAM_INVALID_POS;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_uint_field(Line, &tmp32);
		if (Line != NULL && *Line == '\t') {
			++Line;
			Read->PosQuality = (uint8_t)tmp32;
		} else ret = ERR_SAM_INVALID_MAPQ;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_string_field(Line, &Read->Extension->CIGAR, &tmpStringLen);
		if (Line != NULL && *Line == '\t')
			++Line;
		else ret = ERR_SAM_INVALID_CIGAR;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_string_field(Line, &Read->Extension->RNext, &tmpStringLen);
		if (Line != NULL && *Line == '\t')
			++Line;
		else ret = ERR_SAM_INVALID_RNEXT;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_uint_field(Line, &tmp32);
		if (Line != NULL && *Line == '\t') {
			++Line;
			Read->Extension->PNext = tmp32;
		} else ret = ERR_SAM_INVALID_PNEXT;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_int_field(Line, &tmp32);
		if (Line != NULL && *Line == '\t') {
			++Line;
			Read->Extension->TLen = tmp32;
		} else ret = ERR_SAM_INVALID_TLEN;
	}

	if (ret == ERR_SUCCESS) {
		size_t len = 0;

		Line = _sam_read_string_field(Line, &Read->ReadSequence, &len);
		Read->ReadSequenceLen = (uint32_t)len;
		if (Line != NULL && *Line == '\t') {
			++Line;
		} else ret = ERR_SAM_INVALID_SEQ;
	}

	if (ret == ERR_SUCCESS) {
		Line = _sam_read_string_field(Line, (char **)&Read->Quality, &tmpStringLen);
		if (Line == NULL)
			ret = ERR_SAM_INVALID_QUAL;
	}

	if (ret == ERR_SUCCESS) {
		if (Read->ReadSequenceLen == tmpStringLen)
			read_quality_decode(Read);
		else ret = ERR_SAM_SEQ_QUAL_LEN_MISMATCH;
	}

	if (ret != ERR_SUCCESS) {
		if (Read->Quality != NULL)
			utils_free(Read->Quality);

		if (Read->ReadSequence != NULL)
			utils_free(Read->ReadSequence);

		if (Read->Extension->RNext != NULL)
			utils_free(Read->Extension->RNext);
		
		if (Read->Extension->CIGAR != NULL)
			utils_free(Read->Extension->CIGAR);
			
		if (Read->Extension->RName != NULL)
			utils_free(Read->Extension->RName);

		if (Read->Extension->TemplateName != NULL)
			utils_free(Read->Extension->TemplateName);
	}

	return ret;
}

ERR_VALUE read_create_from_fasta_seq(const char *Seq, const size_t SeqLen, const char *SeqName, const size_t SeqNameLen, PONE_READ *Read)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	return ret;
}


void read_destroy(PONE_READ Read)
{
	_read_destroy_structure(Read);
	utils_free(Read);

	return;
}


void read_set_stats(const ONE_READ *Reads, const size_t Count, const uint8_t MinPosQuality, PBAD_READS_STATISTICS Stats)
{
	memset(Stats, 0, sizeof(BAD_READS_STATISTICS));
	Stats->Total = Count;
	for (size_t i = 0; i < Count; ++i) {
		if (Reads->Extension->Flags.Bits.Paired)
			++Stats->Paired;
		
		if (Reads->Pos == UINT64_MAX ||
			Reads->PosQuality < MinPosQuality ||
			Reads->Extension->Flags.Bits.Unmapped ||
			Reads->Extension->Flags.Bits.Supplementary ||
			Reads->Extension->Flags.Bits.SecondaryAlignment ||
			Reads->Extension->Flags.Bits.Duplicate) {
			++Stats->BadTotal;
			if (Reads->Pos == (uint64_t)-1LL)
				++Stats->BadPosZero;

			if (Reads->PosQuality < MinPosQuality)
				++Stats->BadPosQuality;

			if (Reads->Extension->Flags.Bits.Duplicate)
				++Stats->BadDuplicate;

			if (Reads->Extension->Flags.Bits.Supplementary)
				++Stats->BadSupplementary;

			if (Reads->Extension->Flags.Bits.SecondaryAlignment)
				++Stats->BadSecondaryAlignment;

			if (Reads->Extension->Flags.Bits.Unmapped)
				++Stats->BadUnmapped;
		} else {
			boolean softClipped = FALSE;
			boolean hardClipped = FALSE;
			const char *cigar = Reads->Extension->CIGAR;
			const size_t cigarLen = strlen(cigar);
			
			for (size_t j = 0; j < cigarLen; ++j) {
				switch (cigar[j]) {
					case 'H':
						hardClipped = TRUE;
						break;
					case 'S':
						softClipped = TRUE;
						break;
				}
			}

			if (hardClipped && softClipped)
				++Stats->BothClippedGood;
			else if (hardClipped)
				++Stats->HardClippedGood;
			else if (softClipped)
				++Stats->SoftClippedGood;
		}

		++Reads;
	}


	return;
}


void read_set_stats_print(FILE *Stream, const BAD_READS_STATISTICS *Stats)
{
	fprintf(Stream, "Total reads:         %zu\n", Stats->Total);
	fprintf(Stream, "Total paired:        %zu\n", Stats->Paired);
	fprintf(Stream, "Bad reads:           %zu (%.2lf %%)\n", Stats->BadTotal, (double)Stats->BadTotal * 100 / Stats->Total);
	fprintf(Stream, "Zero POS:            %zu (%.2lf %%)\n", Stats->BadPosZero, (double)Stats->BadPosZero * 100 / Stats->BadTotal);
	fprintf(Stream, "Bad MAPQ:            %zu (%.2lf %%)\n", Stats->BadPosQuality, (double)Stats->BadPosQuality * 100 / Stats->BadTotal);
	fprintf(Stream, "Unmapped:            %zu (%.2lf %%)\n", Stats->BadUnmapped, (double)Stats->BadUnmapped * 100 / Stats->BadTotal);
	fprintf(Stream, "Supplementary:       %zu (%.2lf %%)\n", Stats->BadSupplementary, (double)Stats->BadSupplementary * 100 / Stats->BadTotal);
	fprintf(Stream, "Secondary:           %zu (%.2lf %%)\n", Stats->BadSecondaryAlignment, (double)Stats->BadSecondaryAlignment * 100 / Stats->BadTotal);
	fprintf(Stream, "Duplicate:           %zu (%.2lf %%)\n", Stats->BadDuplicate, (double)Stats->BadDuplicate * 100 / Stats->BadTotal);
	fprintf(Stream, "Soft clipped:        %zu (%.2lf %%)\n", Stats->SoftClippedGood, (double)Stats->SoftClippedGood * 100 / (Stats->Total - Stats->BadTotal));
	fprintf(Stream, "Hard clipped:        %zu (%.2lf %%)\n", Stats->HardClippedGood, (double)Stats->HardClippedGood * 100 / (Stats->Total - Stats->BadTotal));
	fprintf(Stream, "Both clipped:        %zu (%.2lf %%)\n", Stats->BothClippedGood, (double)Stats->BothClippedGood * 100 / (Stats->Total - Stats->BadTotal));

	return;
}


void read_set_destroy(PONE_READ ReadSet, const size_t Count)
{
	PONE_READ tmp = ReadSet;

	for (size_t i = 0; i < Count; ++i) {
		_read_destroy_structure(tmp);
		++tmp;
	}

	utils_free(ReadSet);

	return;
}


ERR_VALUE read_set_merge(PONE_READ *Target, const size_t TargetCount, struct _ONE_READ *Source, const size_t SourceCount)
{
	PONE_READ tmp = NULL;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc_ONE_READ(TargetCount + SourceCount, &tmp);
	if (ret == ERR_SUCCESS) {
		memcpy(tmp, *Target, TargetCount*sizeof(ONE_READ));
		memcpy(tmp + TargetCount, Source, SourceCount*sizeof(ONE_READ));		
		utils_free(Source);
		utils_free(*Target);
		*Target = tmp;
	}

	return ret;
}


void read_adjust(PONE_READ Read, const uint64_t RegionStart, const size_t RegionLength)
{
	uint32_t endStripped = 0;
	uint32_t startStripped = 0;

	if (Read->Pos < RegionStart)
		startStripped = (uint32_t)(RegionStart - Read->Pos);

	if (Read->Pos + Read->ReadSequenceLen >= RegionStart + RegionLength)
		endStripped = (uint32_t)(Read->Pos + Read->ReadSequenceLen - RegionStart - RegionLength);

	if (startStripped > 0) {
		Read->Offset = startStripped;
		if (Read->ReadSequenceLen >= startStripped) {
			Read->ReadSequenceLen -= startStripped;
			Read->ReadSequence += startStripped;
			Read->Quality += startStripped;
			Read->Pos += startStripped;
		} else {
			startStripped -= Read->ReadSequenceLen;
			Read->ReadSequenceLen = 0;
		}
	}

	if (endStripped > 0) {
		if (Read->ReadSequenceLen > endStripped) {
			Read->ReadSequenceLen -= endStripped;
		} else {
			endStripped -= Read->ReadSequenceLen;
			Read->ReadSequenceLen = 0;
		}
	}

	return;
}


void read_split(PONE_READ Read)
{
	READ_PART part;
	boolean end = FALSE;

	if (Read->Extension->CIGAR != NULL && *Read->Extension->CIGAR != '\0' && *Read->Extension->CIGAR != '*') {
		char t;
		unsigned long count = 1;
		const char *c = Read->Extension->CIGAR;

		part.ReadSequence = Read->ReadSequence;
		part.Position = Read->Pos;
		part.Offset = 0;
		part.ReadSequenceLength = 0;
		part.Quality = Read->Quality;
		while (*c != '\0') {
			char *tmp;

			assert(!end);
			count = strtoul(c, &tmp, 10);
			c = tmp;
			t = *c;
			if (t != '\0' && count > 0) {
				++c;
				switch (t) {
				case 'M':
				case 'I':
					part.ReadSequenceLength += count;
					break;
				case 'D':
					break;
				case 'S':
					if (part.ReadSequenceLength > 0) {
						end = TRUE;
					} else {
						part.ReadSequence += count;
						part.Quality += count;
					}
					break;
				case 'H':
					if (part.ReadSequenceLength > 0)
						end = TRUE;
					break;
				default:
					part.Offset = 0;
					part.Position = Read->Pos;
					part.ReadSequence = Read->ReadSequence;
					part.ReadSequenceLength = Read->ReadSequenceLen;
					part.Quality = Read->Quality;
					break;
				}
			}
		}

		if (part.ReadSequenceLength > 0 && part.ReadSequenceLength != Read->ReadSequenceLen) {
			Read->Pos = part.Position;
			memmove(Read->ReadSequence, part.ReadSequence, part.ReadSequenceLength*sizeof(char));
			memmove(Read->Quality, part.Quality, part.ReadSequenceLength*sizeof(char));
			Read->ReadSequenceLen = part.ReadSequenceLength;
			Read->ReadSequence[Read->ReadSequenceLen] = '\0';
			Read->Quality[Read->ReadSequenceLen] = '\0';
		}

		if (Read->Extension->CIGAR != NULL) {
			Read->Extension->CIGAR[0] = '*';
			Read->Extension->CIGAR[1] = '\0';
		}
	}

	return;
}


void read_shorten(PONE_READ Read, const uint32_t Count)
{
	if (Read->ReadSequenceLen > 2 * Count) {
		if (!Read->NoEndStrip)
			Read->ReadSequenceLen -= Count;

		Read->ReadSequence[Read->ReadSequenceLen] = '\0';
		Read->Quality[Read->ReadSequenceLen] = '\0';
	}

	return;
}


ERR_VALUE read_append(PONE_READ Read, const char *Seq, const uint8_t *Quality, size_t Length)
{
	char *tmpSeq = NULL;
	uint8_t *tmpQ = NULL;
	const size_t totaLLength = Read->ReadSequenceLen + Length;
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_calloc_char(totaLLength + 1, &tmpSeq);
	if (ret == ERR_SUCCESS) {
		ret = utils_calloc_uint8_t(totaLLength + 1, &tmpQ);
		if (ret == ERR_SUCCESS) {
			memcpy(tmpSeq, Read->ReadSequence, Read->ReadSequenceLen * sizeof(char));
			memcpy(tmpSeq + Read->ReadSequenceLen, Seq, Length * sizeof(char));
			tmpSeq[totaLLength] = '\0';
			memcpy(tmpQ, Read->Quality, Read->ReadSequenceLen * sizeof(char));
			memcpy(tmpQ + Read->ReadSequenceLen, Quality, Length * sizeof(uint8_t));
			tmpQ[totaLLength] = '\0';
			utils_free(Read->ReadSequence);
			Read->ReadSequence = tmpSeq;
			utils_free(Read->Quality);
			Read->Quality = tmpQ;
			Read->ReadSequenceLen = totaLLength;
			if (ret != ERR_SUCCESS)
				utils_free(tmpQ);
		}

		if (ret != ERR_SUCCESS)
			utils_free(tmpSeq);
	}

	return ret;
}



/************************************************************************/
/*                ASSEMBLY TATKS                                        */
/************************************************************************/

void assembly_task_init(PASSEMBLY_TASK Task, const char *RefSeq, const size_t RefSeqLen, const char *Alternate1, const size_t Alternate1Length, const char *Alternate2, const size_t Alternate2Length, const ONE_READ *ReadSet, const size_t ReadCount)
{
	memset(Task, 0, sizeof(ASSEMBLY_TASK));
	Task->Allocated = FALSE;
	Task->Reference = RefSeq;
	Task->ReferenceLength = RefSeqLen;
	Task->Alternate1 = Alternate1;
	Task->Alternate2 = Alternate2;
	Task->Alternate1Length = Alternate1Length;
	Task->Alternate2Length = Alternate2Length;
	Task->Reads = ReadSet;
	Task->ReadCount = ReadCount;
	
	return;
}


void assembly_task_finit(PASSEMBLY_TASK Task)
{
	if (Task->Allocated) {
		read_set_destroy((PONE_READ)Task->Reads, Task->ReadCount);
		utils_free((char *)Task->Alternate2);
		utils_free((char *)Task->Alternate1);
		utils_free((char *)Task->Reference);
	}

	return;
}


void assembly_task_set_name(PASSEMBLY_TASK Task, const char *Name)
{
	Task->Name = Name;

	return;
}
