#ifndef DEBUG_INFO_H
#define DEBUG_INFO_H

#include <stdio.h>
#include "cfg_decoding.h"

#define CFG_DEBUG_SYS			1
#if (1 == TEST_MODE)
#define CFG_DEBUG_IMPOTANT		1
#define CFG_DEBUG_INFO			1
#define CFG_DEBUG_NOTIE			1
#else
#define CFG_DEBUG_IMPOTANT		0
#define CFG_DEBUG_INFO			0
#define CFG_DEBUG_NOTICE		0
#endif

#if (0 == DEBUG_LOG)

#if (1 == CFG_DEBUG_SYS)
#define DEBUG_SYS				printf
#else
#define DEBUG_SYS				//printf
#endif
#if (1 == CFG_DEBUG_IMPOTANT)
#define DEBUG_IMPOTANT			printf
#else
#define DEBUG_IMPOTANT			//printf
#endif
#if (1 == CFG_DEBUG_INFO)
#define DEBUG_INFO				printf
#else
#define DEBUG_INFO				//printf
#endif
#if (1 == CFG_DEBUG_NOTICE)
#define DEBUG_NOTICE			printf
#else
#define DEBUG_NOTICE			//printf
#endif

#else

#define DEBUG_SYS			fprintf(frc_debug, 
#define DEBUG_IMPOTANT		fprintf(frc_debug, 
#define DEBUG_INFO			fprintf(frc_debug, 
#define DEBUG_NOTICE		fprintf(frc_debug, 

#endif

#endif
