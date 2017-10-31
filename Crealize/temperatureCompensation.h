//*****************************************************************************
//		Copyright (C), 1993-2013, zhejiang supcon instrument Co., Ltd.
//
//		文 件 名：			TemperatureCompensation.h
//
//		功    能：			温度补偿头文件
//
//		版    本：			V1.00
//		作    者：		    吴静超
//		创建日期：			2017.10.26
//
//		修改次数：
//		修改时间：
//		修改内容：
//		修 改 者：
//*****************************************************************************
/**
* 包含库头文件
*/
#include <stdint.h>

//补偿类型分类
#define CRUDE_OIL_TYPE      0   //原油
#define PRODUCT_OIL_TYPE    1   //产品油
#define LUBE_OIL_TYPE       2   //润滑油

/**
*   对外开放接口
*/
//计算补偿后的密度值
uint8_t CalCompenDensity( uint8_t oilKind, double refDensity, double refTemper, double realTemper, double * pRetDensity, double * pVCF20 );

//test
//double getRealRou15(uint8_t oilType, float tstDen, float tstTemp );
