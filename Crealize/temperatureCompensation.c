//*****************************************************************************
//		Copyright (C), 1993-2013, zhejiang supcon instrument Co., Ltd.
//
//		文 件 名：			TemperatureCompensation.c
//
//		功    能：			温度补偿文件
//      说    明：          密度单位为g/m3或kg/L
//
//		版    本：			V1.00�b
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
#include <stdio.h>
#include <stdint.h>
#include <math.h>

/**
* 包含内部头文件
*/
#include "temperatureCompensation.h"

/**
* 内部宏定义
*/
#define DBG_PRINT  printf

//查表区域区分
#define     RAW_OIL          (0x10)     //原油
#define     PETROL           (0x11)     //汽油
#define     UNDUE_OIL        (0x12)     //过渡区
#define     SPACE_OIL        (0x13)     //航煤
#define     DIESEL_OIL       (0x14)     //柴油
#define     LUBE_OIL         (0x15)     //润滑油

#define     ITER_PRECISION         (0.000001)  //VCF计算时的迭代精度
#define     ITER_MAX                99        //迭代最大次数


/******************************************************************************************
 * @fn           static double getBuChang(float realTemp, float stdTemp)
 * @Description  获得补偿系数
 * @Parameters   float realTemp：实际温度
                 float stdTemp：标准温度
 * @Rerurn
 * @Notes
 * @Author       WYB
 * @Creatdate    2016-09-22
 *******************************************************************************************/
static double getComModulus(float realTemp, float stdTemp)
{
	double ret;
	double rt, st;

	rt = realTemp;
	st = stdTemp;

	ret = 1-2.3*pow(10,-5)*(rt-st)-2*pow(10,-8)*pow((rt-st),2);
    return ret;
}

/******************************************************************************************
 * @fn           static double getAlpha(uint8_t oilType, double density)
 * @Description  计算阿尔法系数
 * @Parameters   uint8_t oilType ：油品类型
                 double density：密度
 * @Rerurn
 * @Notes
 * @Author       WYB
 * @Creatdate    2016-09-22
 *******************************************************************************************/
static double getAlpha(uint8_t oilType, double density)
{
	double K0,K1,A;
	double ret;

	switch(oilType)
	{
		case RAW_OIL:  //原油
				K0 = 613.9723;
				K1 = 0;
				A = 0;
			break;
		case PETROL:   //汽油
				K0 = 346.4228;
				K1 = 0.4388;
				A = 0;
			break;
		case UNDUE_OIL:   //过渡区
				K0 = 2680.321;
				K1 = 0;
				A = -0.00336312;
			break;
		case SPACE_OIL:   //航煤
				K0 = 594.5418;
				K1 = 0;
				A = 0;
			break;
		case DIESEL_OIL:  //柴油
				K0 = 186.9696;
				K1 = 0.4862;
				A = 0;
			break;
		case LUBE_OIL:    //润滑油
				K0 = 0;
				K1 = 0.6278;
				A = 0;
			break;
		default:
				K0 = 346.4228;
				K1 = 0.4388;
				A = 0;
			break;
	}

	ret = K0/(density*density) + K1/density + A;
	return ret;

}

/******************************************************************************************
 * @fn           static double getExp(double alpha,double rou,double deltTemp)
 * @Description
 * @Parameters

 * @Rerurn
 * @Notes
 * @Author       WYB
 * @Creatdate    2016-09-22
 *******************************************************************************************/
static double getExp(double alpha,double deltTemp)
{
	double ret;

	ret = exp(-1*alpha*deltTemp*(1+0.8*alpha*deltTemp));
	return ret;
}

/** \brief
 *
 * \param oil uint8_t
 * \return uint8_t
 *
 */
static uint8_t convToCompenChart( uint8_t oil )
{
    uint8_t ret;

    switch( oil )
    {
        case CRUDE_OIL_TYPE:
        {
            ret = RAW_OIL;
            break;
        }
        case PRODUCT_OIL_TYPE:
        {
            ret = PETROL;       //若使用成品油，先计算汽油的15℃标准密度，再后续再选择表格
            break;
        }
        case LUBE_OIL_TYPE:
        {
            ret = LUBE_OIL;
            break;
        }
        default:
        {
            ret = 0;
        }
    }

    return ret;
}

/** \brief 根据汽油的15℃标准密度确定补偿表格
 *
 * \param den15 double
 * \return uint8_t
 *
 */
static uint8_t getChartForProductOil( double den15 )
{
    uint8_t ret = 0;

    if( (den15>=838.3) && (den15<1163.5) )          //柴油
    {
        ret = DIESEL_OIL;
    }
    else if( (den15>=787.5) && (den15<838.3) )      //航煤
    {
        ret = SPACE_OIL;
    }
    else if( (den15>=770.3) && (den15<787.5) )      //过渡区
    {
        ret = UNDUE_OIL;
    }
    else if( (den15>=610.6) && (den15<770.3) )      //汽油
    {
        ret = PETROL;
    }

    return ret;

}

/******************************************************************************************
 * @fn           static double getRealRou15(uint_8 oilType, double tstDen, double tstTemp )
 * @Description  迭代法求解15℃时的标准密度
 * @Parameters   uint_8 oilType：油品类型
                 double tstDen：实验密度，视密度
                 double tstTemp：实验温度
 * @Rerurn
 * @Notes
 * @Author       WYB
 * @Creatdate    2016-09-22
 *******************************************************************************************/
static double getRealRou15(uint8_t oilType, double tstDen, double tstTemp )
{
	double rou15,tempRou;
	double alpha;
	double rou_t;

	uint16_t count=0;

	tempRou = tstDen;

    do{
        rou15 = tempRou;
        rou_t = tstDen*getComModulus(tstTemp, 15);
        alpha = getAlpha(oilType, rou15);
        tempRou = rou_t/getExp(alpha,(double)(tstTemp-15));
        count++;
        if(count%10==0)
        {
            #ifdef EN_WDOG
                RstWdt();
            #endif
        }
        else if(count > ITER_MAX)
        {
            DBG_PRINT("-->getRealRou15:Iterate too much! \n");
            break;
        }
        DBG_PRINT("-->getRealRou15:rou15=%lf, tempRou=%lf, alpha=%lf!\n",rou15,tempRou,alpha);
	}while(fabs(rou15-tempRou)>ITER_PRECISION);
	DBG_PRINT("-->getRealRou15:count=%3d, rou15=%f!\n",count,tempRou);
	return tempRou;
}

/******************************************************************************************
 * @fn           static double getRealRou20(uint8_t oilType, double rou15)
 * @Description  获取20℃油品密度
 * @Parameters   uint8_t oilType：油品类型
                 double rou15：油品15℃的标准密度

 * @Rerurn       20℃标准密度，单精度浮点数据
 * @Notes
 * @Author       WYB
 * @Creatdate    2016-09-22
 *******************************************************************************************/
static double getRealRou20(uint8_t oilType, double rou15)
{
	double rou20;
	double alpha;
	double VCF15;

	alpha = getAlpha(oilType, rou15);       //TODO此处为重复计算
	VCF15 = getExp(alpha,(20-15));
	rou20 = rou15 * VCF15;
	DBG_PRINT("-->getRealRou20:Raw rou20=%f!\n",rou20);
	return rou20;
}

/** \brief
 *
 * \param oilType uint8_t
 * \param realTemper double
 * \param den15 double
 * \param den20 double
 * \return double
 *
 */
static double getVCF20( uint8_t type, double temper, double den15, double den20 )
{
    double ret, alpha;

    alpha = getAlpha(type, den15);    //TODO此处为重复计算
    ret = den15*getExp(alpha,(temper-15));
    ret = ret/den20;

    return ret;
}


/** \brief
 *
 * \param oilKind uint8_t
 * \param refDensity double
 * \param refTemper double
 * \param realTemper double
 * \param pRetDensity double*
 * \param pVCF20 double*
 * \return uint8_t
 *
 */
uint8_t CalCompenDensity( uint8_t oilKind, double refDensity, double refTemper, double realTemper, double * pRetDensity, double * pVCF20 )
{
    double den15;
    double retDen, retVCF20;
    uint8_t ret = 0, compenChart;

    //对油品进行检查
    if( (oilKind < CRUDE_OIL_TYPE) || (oilKind > LUBE_OIL_TYPE) )
    {
        return ret;
    }
    compenChart = convToCompenChart( oilKind );
    if( 0 == compenChart )  //油类型有误
    {
        return ret;
    }

    if( PETROL == compenChart )
    {
        den15 = getRealRou15( compenChart, refDensity, refTemper );
        compenChart = getChartForProductOil( den15 );
        if( 0 == compenChart )  //成品油密度空间有误
        {
            return ret;
        }
    }

    den15 = getRealRou15( compenChart, refDensity, refTemper );
    retDen = getRealRou20( compenChart, den15 );
    retVCF20 = getVCF20( compenChart, realTemper, den15, retDen );

    *pRetDensity = retDen*getComModulus( 15, 20 );  //计算20度标准密度真值;
    *pVCF20 = retVCF20;
    ret = 1;

    return ret;
}

