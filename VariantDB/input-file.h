
#ifndef __GASSM_INPUT_FILE_H__
#define __GASSM_INPUT_FILE_H__


#include "gen_dym_array.h"
#include "pointer_array.h"
#include "reads.h"


typedef enum _EActiveRegionType {
	artUnknown,
	artValid,
} EActiveRegionType, *PEActiveRegionType;

typedef struct _ACTIVE_REGION {
	EActiveRegionType Type;
	uint64_t Offset;
	uint64_t Length;
	const char *Sequence;
} ACTIVE_REGION, *PACTIVE_REGION;

GEN_ARRAY_TYPEDEF(ACTIVE_REGION);
GEN_ARRAY_IMPLEMENTATION(ACTIVE_REGION)
POINTER_ARRAY_TYPEDEF(ACTIVE_REGION);
POINTER_ARRAY_IMPLEMENTATION(ACTIVE_REGION)

typedef struct _FASTA_FILE {
	char *FileData;
	size_t DataLength;
	char *CurrentPointer;
} FASTA_FILE, *PFASTA_FILE;

typedef struct _REFSEQ_DATA {
	const char *Sequence;
	size_t Length;
	uint64_t StartPos;
	const char *Name;
} REFSEQ_DATA, *PREFSEQ_DATA;

typedef ERR_VALUE (INPUT_READ_CALLBACK)(const ONE_READ *Read, void *Context);

typedef enum _EVCFVariantType {
	vcfvtUnknown,
	vcfvtSNP,
	vcfvtInsertion,
	vcfvtDeletion,
	vcfvtReplace,
	vcfvtMax,
} EVCFVariantType, *PEVCFVariantType;

typedef struct _VCF_VARIANT {
	char *Chrom;
	unsigned long long Pos;
	char *ID;
	char *Ref;
	char *Alt;
	unsigned long Quality;
	size_t ReadSupport;
	size_t TotalReadsAtPosition;
	EVCFVariantType Type;
	struct _VCF_VARIANT *Alternative;
} VCF_VARIANT, *PVCF_VARIANT;


GEN_ARRAY_TYPEDEF(VCF_VARIANT);
GEN_ARRAY_IMPLEMENTATION(VCF_VARIANT)
POINTER_ARRAY_TYPEDEF(VCF_VARIANT);
POINTER_ARRAY_IMPLEMENTATION(VCF_VARIANT)

typedef struct _CONFIDENT_REGION {
	char *Chrom;
	unsigned long long Start;
	unsigned long long End;
} CONFIDENT_REGION, *PCONFIDENT_REGION;

GEN_ARRAY_TYPEDEF(CONFIDENT_REGION);
GEN_ARRAY_IMPLEMENTATION(CONFIDENT_REGION)
POINTER_ARRAY_TYPEDEF(CONFIDENT_REGION);
POINTER_ARRAY_IMPLEMENTATION(CONFIDENT_REGION)

typedef struct _VCF_VARIANT_FILDER {
	const CONFIDENT_REGION *Regions;
	size_t RegionCount;
} VCF_VARIANT_FILTER, *PVCF_VARIANT_FILTER;


ERR_VALUE fasta_load(const char *FileName, PFASTA_FILE FastaRecord);
ERR_VALUE fasta_read_seq(PFASTA_FILE FastaRecord, PREFSEQ_DATA Data);
void fasta_free_seq(PREFSEQ_DATA Data);
void fasta_free(PFASTA_FILE FastaRecord);

ERR_VALUE input_get_reads(const char *Filename, const CONFIDENT_REGION *Region, INPUT_READ_CALLBACK *Callback, void *Context);

ERR_VALUE input_refseq_to_regions(const char *RefSeq, const size_t RefSeqLen, PACTIVE_REGION *Regions, size_t *Count);
ERR_VALUE input_get_region_by_offset(const PACTIVE_REGION Regions, const size_t Count, const uint64_t Offset, size_t *Index, uint64_t *RegionOffset);
void input_free_regions(PACTIVE_REGION Regions, const size_t Count);

ERR_VALUE input_variant_create(const char *Chrom, const char *ID, unsigned long long Pos, const char *Ref, const char *Alt, unsigned long Quality, PVCF_VARIANT Variant);
ERR_VALUE input_get_variants(const char *FileName, const VCF_VARIANT_FILTER *Filter, PGEN_ARRAY_VCF_VARIANT Array);
void input_free_variant(const VCF_VARIANT *Variant);
void input_Free_variants(PGEN_ARRAY_VCF_VARIANT Array);
boolean input_variant_in_filter(const VCF_VARIANT_FILTER *Filter, const VCF_VARIANT *Variant);
boolean input_variant_normalize(const char *Reference, PVCF_VARIANT Variant);
boolean input_variant_equal(const VCF_VARIANT *A, const VCF_VARIANT *B);

ERR_VALUE input_get_bed(const char *FileName, const CONFIDENT_REGION *Area, PGEN_ARRAY_CONFIDENT_REGION Array);
void input_free_bed(PGEN_ARRAY_CONFIDENT_REGION Array);



#endif 
