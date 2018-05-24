
#ifndef __VARIANT_DB_H__
#define __VARIANT_DB_H__


#define VDB_OPTION_REF_FILE				"ref-file"
#define VDB_OPTION_SAM_FILE				"sam-file"
#define VDB_OPTION_VCF_FILE				"vcf-file"
#define VDB_OPTION_BED_FILE				"bed-file"
#define VDB_OPTION_CHROM				"chrom"
#define VDB_OPTION_START				"start"
#define VDB_OPTION_STOP					"stop"
#define VDB_OPTION_DONT_NORMALIZE		"dont-normalize"
#define VDB_OPTION_MAX_MS				"max-ms"
#define VDB_OPTION_HELP					"help"
#define VDB_OPTION_VERBOSE				"verbose"

#define VDB_OPTION_REF_FILE_DESC		"ref-file"
#define VDB_OPTION_SAM_FILE_DESC		"sam-file"
#define VDB_OPTION_VCF_FILE_DESC		"vcf-file"
#define VDB_OPTION_BED_FILE_DESC		"bed-file"
#define VDB_OPTION_CHROM_DESC			"chrom"
#define VDB_OPTION_START_DESC			"start"
#define VDB_OPTION_STOP_DESC			"stop"
#define VDB_OPTION_DONT_NORMALIZE_DESC	"dont-normalize"
#define VDB_OPTION_MAX_MS_DESC			"max-ms"
#define VDB_OPTION_HELP_DESC			"help"
#define VDB_OPTION_VERBOSE_DESC			"verbose"

#define VDB_OPTION_REF_FILE_SHORT		'f'
#define VDB_OPTION_SAM_FILE_SHORT		's'
#define VDB_OPTION_VCF_FILE_SHORT		'v'
#define VDB_OPTION_BED_FILE_SHORT		'b'
#define VDB_OPTION_CHROM_SHORT			'c'
#define VDB_OPTION_START_SHORT			'f'
#define VDB_OPTION_STOP_SHORT			't'
#define VDB_OPTION_HELP_SHORT			'h'
#define VDB_OPTION_VERBOSE_SHORT		'V'
#define VDB_OPTION_DONT_NORMALIZE_SHORT	'n'
#define VDB_OPTION_MAX_MS_SHORT			'm'



#define CMD_OPTION_INIT(aMacroName, aType, aDefault)	\
	option_add_##aType(aMacroName, aDefault);	\
	option_set_shortcut(aMacroName, aMacroName##_SHORT);	\
	option_set_description_const(aMacroName, aMacroName##_DESC);	\

#define CMD_OPTION_GET(aMacroName, aType, aPResult)	\
	option_get_##aType(aMacroName, aPResult);	\

#endif
