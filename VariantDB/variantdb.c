

#include <stdio.h>
#include "err.h"
#include "utils.h"
#include "options.h"




int main(int argc, char **argv)
{
	ERR_VALUE ret = ERR_INTERNAL_ERROR;

	ret = utils_allocator_init(1);
	if (ret == ERR_SUCCESS) {
		ret = options_module_init(37);
		if (ret == ERR_SUCCESS) {

			options_module_finit();
		}
	}

	return 0;
}
