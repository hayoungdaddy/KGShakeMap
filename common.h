#ifndef COMMON_H
#define COMMON_H

#include <QString>
#include <QDateTime>
#include <QFile>
#include <QTextStream>
#include <QDir>
#include <QDebug>
#include <QColor>
#include <QGradient>
#include <QtMath>

#include <proj_api.h>
#include <math.h>

#define XWIDTH  834896.18
#define YHEIGHT 656784.99
#define XORIGIN 13859276.60
#define YORIGIN 3673543.20

#define SCREEN_WIDTH    814
#define SCREEN_HEIGHT   768
#define SCREEN_BPP      32

#define MAX_DISTANCE    10000

#define EQ_DEPTH 15
#define VALID_RANGE 0.7
#define VALID_RANGE2 0.7

typedef struct _xyland
{
    QList<int> xList;
    QList<int> yList;
} XYLAND;

typedef struct _event
{
    QString evid;
    QDateTime origintime;
    double lat;
    double lon;
    double depth;
    double ml;
    double mb;
    QString loc;
    int nSta;
} _EVENT;

typedef struct _station
{
    int index;
    QString staName;
    double lat;
    double lon;
    double distance;
    int mapX;
    int mapY;
    float pga;
    double predPGA;
} _STATION;

typedef struct _point
{
    int index;
    int landX;
    int landY;
    double mapZ;
    QList<_STATION> staList;
    QList<double> mapLUT;
    float predPGA;
    float realPGA;
} _POINT;

static double getPredictedValue(double dist, double mag, double depth)
{
    double y = (dist * dist) + (depth * depth);
    double R = sqrt(y);
    double c0, c1, c2;

    /* 2001  */
    double ksai0[4] = {0.1250737E+02, 0.4874629E+00, -0.2940726E-01, 0.1737204E-01};
    double ksai1[4] = {-0.1928185E-02, 0.2251016E-03, -0.6378615E-04, 0.6967121E-04};
    double ksai2[4] = {-0.5795112E+00, 0.1138817E+00, -0.1162326E-01, -0.3646674E-02};

    c0 = ksai0[0] + ksai0[1]*(mag-6) + ksai0[2]*pow((mag-6), 2) + ksai0[3]*pow((mag-6), 3);
    c1 = ksai1[0] + ksai1[1]*(mag-6) + ksai1[2]*pow((mag-6), 2) + ksai1[3]*pow((mag-6), 3);
    c2 = ksai2[0] + ksai2[1]*(mag-6) + ksai2[2]*pow((mag-6), 2) + ksai2[3]*pow((mag-6), 3);

    /* 2003
    double ksai0[4] = {0.1073829E+02 , 0.5909022E+00, -0.5622945E-01, 0.2135007E-01};
    double ksai1[4] = {-0.2379955E-02, 0.2081359E-03, -0.2046806E-04, 0.4192630E-04};
    double ksai2[4] = {-0.2437218E+00, 0.9498274E-01, -0.8804236E-02, -0.3302350E-02};

    c0 = ksai0[0] + (ksai0[1]*(mag-6)) + pow((ksai0[2]*(mag-6)), 2) + pow((ksai0[3]*(mag-6)), 3);
    c1 = ksai1[0] + (ksai1[1]*(mag-6)) + pow((ksai1[2]*(mag-6)), 2) + pow((ksai1[3]*(mag-6)), 3);
    c2 = ksai2[0] + (ksai2[1]*(mag-6)) + pow((ksai2[2]*(mag-6)), 2) + pow((ksai2[3]*(mag-6)), 3);
     */

    double lnSA;
    if(R < 100)
        lnSA = c0 + c1*R + c2*log(R) - log(R) - 0.5*log(100);
    else if(R > 100)
        lnSA = c0 + c1*R + c2*log(R) - log(100) - 0.5*log(R);
    else if(R == 100)
        lnSA = c0 + c1*R + c2*log(R) - log(100) - 0.5*log(100);

    //qDebug() << c0 << c1 << c2 << R << exp(lnSA) << qLn(lnSA);

    return exp(lnSA);
}

static void getMapZValue(_POINT &point)
{
    double A=0., B=0., z=0., W=0.;

    for(int i=0;i<point.staList.size();i++)
    {
        if(point.landX == point.staList.at(i).mapX && point.landY == point.staList.at(i).mapY)
        {
            A = point.staList.at(i).pga;
            B = 1;
            break;
        }

        z = point.staList.at(i).pga;
        W = point.mapLUT.at(i);

        A += ( W * z );
        B += ( W );

    }
    point.realPGA = A/B;

    if(point.realPGA != point.realPGA)
        point.realPGA = 0;

    /*
    if(point.predPGA - (point.predPGA * VALID_RANGE2) <= point.realPGA &&
        point.predPGA + (point.predPGA * VALID_RANGE2) >= point.realPGA)
        point.mapZ = point.realPGA;
    else
        point.mapZ = point.predPGA;
        */


    // using prediction value
    //point.mapZ = point.predPGA;

    // using real value
    point.mapZ = point.realPGA;
};

static int redColor(float gal)
{
    int color ;

    // red color value
    if( gal <= 0.0098 )
    {
      color = 191 ;
    }
    else if (gal > 0.0098 && gal <= 0.0392)
    {
      color = gal * (-3265.31) + 223 ;
    }
    else if (gal > 0.0392 && gal <= 0.0784)
    {
      color = 95 ;
    }
    else if (gal > 0.0784 && gal <= 0.098)
    {
      color = gal * 3265.31 - 161 ;
    }
    else if (gal > 0.098 && gal <= 0.98)
    {
      color = gal * 103.82 + 148.497 ;
    }
    else if (gal > 0.98 && gal <= 147)
    {
      color = 255 ;
    }
    else if (gal > 147 && gal <= 245)
    {
      color = -0.00333195 * pow(gal,2) + 0.816327 * gal + 207 ;
    }
    else if (gal > 245)
      color = 207 ;

    return color ;
}

static int greenColor(float gal)
{
    int color ;
    // red color value
    if( gal <= 0.98 )
    {
      color = 255 ;
    }
    else if (gal > 0.98 && gal <= 9.8)
    {
      color = -0.75726 * gal * gal - 0.627943 * gal + 255.448 ;
    }
    else if (gal > 0.98 && gal <= 245)
    {
      color = 0.00432696 * gal * gal - 1.84309 * gal + 192.784 ;
      if(color < 0)
        color = 0 ;
    }
    else if (gal > 245)
      color = 0 ;

    return color ;
}

static int blueColor(float gal)
{
    int color ;

    // red color value
    if( gal <= 0.0098 )
    {
      color = 255 ;
    }
    else if (gal > 0.0098 && gal <= 0.098)
    {
      color = -19799.2 * gal * gal + 538.854 * gal + 260.429 ;
    }
    else if (gal > 0.098 && gal <= 0.98)
    {
      color = -35.4966 * gal * gal - 65.8163 * gal + 116.264 ;
    }
    else if (gal > 0.98 && gal <= 3.92)
    {
      color = -5.10204 * gal + 20 ;
    }
    else if (gal > 3.92)
    {
      color = 0 ;
    }

    if(color > 255)
      color = 255 ;

    return color ;
}

#define PI 3.14159265358979323846

static int geo_to_km(double lat1,double lon1,double lat2,double lon2,double* dist,double* azm)
{
    double a, b;
    double semi_major=a=6378.160;
    double semi_minor=b=6356.775;
    double torad, todeg;
    double aa, bb, cc, dd, top, bottom, lambda12, az, temp;
    double v1, v2;
    double fl, e, e2, eps, eps0;
    double b0, x2, y2, z2, z1, u1p, u2p, xdist;
    double lat1rad, lat2rad, lon1rad, lon2rad;
    double coslon1, sinlon1, coslon2, sinlon2;
    double coslat1, sinlat1, coslat2, sinlat2;
    double tanlat1, tanlat2, cosazm, sinazm;

    double c0, c2, c4, c6;

    double c00=1.0, c01=0.25, c02=-0.046875, c03=0.01953125;
    double c21=-0.125, c22=0.03125, c23=-0.014648438;
    double c42=-0.00390625, c43=0.0029296875;
    double c63=-0.0003255208;

    if( lat1 == lat2 && lon1 == lon2 ) {
        *azm = 0.0;
        *dist= 0.0;
        return(1);
    }

    torad = PI / 180.0;
    todeg = 1.0 / torad;
    fl = ( a - b ) / a;
    e2 = 2.0*fl - fl*fl;
    e  = sqrt(e2);
    eps = e2 / ( 1.0 - e2);

    temp=lat1;
    if(temp == 0.) temp=1.0e-08;
    lat1rad=torad*temp;
    lon1rad=torad*lon1;

    temp=lat2;
    if(temp == 0.) temp=1.0e-08;
    lat2rad=torad*temp;
    lon2rad=torad*lon2;

    coslon1 = cos(lon1rad);
    sinlon1 = sin(lon1rad);
    coslon2 = cos(lon2rad);
    sinlon2 = sin(lon2rad);
    tanlat1 = tan(lat1rad);
    tanlat2 = tan(lat2rad);
    sinlat1 = sin(lat1rad);
    coslat1 = cos(lat1rad);
    sinlat2 = sin(lat2rad);
    coslat2 = cos(lat2rad);

    v1 = a / sqrt( 1.0 - e2*sinlat1*sinlat1 );
    v2 = a / sqrt( 1.0 - e2*sinlat2*sinlat2 );
    aa = tanlat2 / ((1.0+eps)*tanlat1);
    bb = e2*(v1*coslat1)/(v2*coslat2);
    lambda12 = aa + bb;
    top = sinlon2*coslon1 - coslon2*sinlon1;
    bottom = lambda12*sinlat1-coslon2*coslon1*sinlat1-sinlon2*sinlon1*sinlat1;
    az = atan2(top,bottom)*todeg;
    if( az < 0.0 ) az = 360 + az;
    *azm = az;
    az = az * torad;
    cosazm = cos(az);
    sinazm = sin(az);

    if( lat2rad < 0.0 )
    {
        temp = lat1rad;
        lat1rad = lat2rad;
        lat2rad = temp;
        temp = lon1rad;
        lon1rad = lon2rad;
        lon2rad = temp;

        coslon1 = cos(lon1rad);
        sinlon1 = sin(lon1rad);
        coslon2 = cos(lon2rad);
        sinlon2 = sin(lon2rad);
        tanlat1 = tan(lat1rad);
        tanlat2 = tan(lat2rad);
        sinlat1 = sin(lat1rad);
        coslat1 = cos(lat1rad);
        sinlat2 = sin(lat2rad);
        coslat2 = cos(lat2rad);

        v1 = a / sqrt( 1.0 - e2*sinlat1*sinlat1 );
        v2 = a / sqrt( 1.0 - e2*sinlat2*sinlat2 );

        aa = tanlat2 / ((1.0+eps)*tanlat1);
        bb = e2*(v1*coslat1)/(v2*coslat2);
        lambda12 = aa + bb;

        top = sinlon2*coslon1 - coslon2*sinlon1;
        bottom =lambda12*sinlat1-coslon2*coslon1*sinlat1-
            sinlon2*sinlon1*sinlat1;
        az = atan2(top,bottom);
        cosazm = cos(az);
        sinazm = sin(az);

    }

    eps0 = eps * ( coslat1*coslat1*cosazm*cosazm + sinlat1*sinlat1 );
    b0 = (v1/(1.0+eps0)) * sqrt(1.0+eps*coslat1*coslat1*cosazm*cosazm);

    x2 = v2*coslat2*(coslon2*coslon1+sinlon2*sinlon1);
    y2 = v2*coslat2*(sinlon2*coslon1-coslon2*sinlon1);
    z2 = v2*(1.0-e2)*sinlat2;
    z1 = v1*(1.0-e2)*sinlat1;

    c0 = c00 + c01*eps0 + c02*eps0*eps0 + c03*eps0*eps0*eps0;
    c2 =       c21*eps0 + c22*eps0*eps0 + c23*eps0*eps0*eps0;
    c4 =                  c42*eps0*eps0 + c43*eps0*eps0*eps0;
    c6 =                                  c63*eps0*eps0*eps0;

    bottom = cosazm*sqrt(1.0+eps0);
    u1p = atan2(tanlat1,bottom);

    top = v1*sinlat1+(1.0+eps0)*(z2-z1);
    bottom = (x2*cosazm-y2*sinlat1*sinazm)*sqrt(1.0+eps0);
    u2p = atan2(top,bottom);

    aa = c0*(u2p-u1p);
    bb = c2*(sin(2.0*u2p)-sin(2.0*u1p));
    cc = c4*(sin(4.0*u2p)-sin(4.0*u1p));
    dd = c6*(sin(6.0*u2p)-sin(6.0*u1p));

    xdist = fabs(b0*(aa+bb+cc+dd));
    *dist = xdist;
    return(1);
}

static double getDistance(double lat1, double lon1, double lat2, double lon2)
{
    double dist, azim;
    int rtn = geo_to_km(lat1, lon1, lat2, lon2, &dist, &azim);
    return dist;
}

static QColor getColorfromGal(float value)
{
    //float pg = value / 980 * 100; // convert gal to %g
    float pg = value;

    QList<float> pgaRange;
    //pgaRange << 0 << 0.1 << 0.3 << 0.5 << 2.4 << 6.7 << 13 << 24 << 44 << 83;
    pgaRange << 0 << 0.98 << 2.94 << 4.90 << 23.52 << 65.66 << 127.40 << 235.20 << 431.20 << 813.40;
    QList<QColor> pgaColor;
    pgaColor << QColor("#FFFFFF") << QColor("#A5DDF9") << QColor("#92D050") << QColor("#FFFF00")
             << QColor("#FFC000") << QColor("#FF0000") << QColor("#A32777") << QColor("#632523")
             << QColor("#4C2600") << QColor("#000000");

    int index;

    for(int i=1;i<9;i++)
    {
        if(pg < 0.98)
        {
            index = 0;
            break;
        }
        else if(pg > 813.40)
        {
            index = 9;
            return pgaColor.at(index);
            break;
        }
        else
        {
            if(pg >= pgaRange.at(i) && pg < pgaRange.at(i+1))
            {
                index = i;
                break;
            }
        }
    }

    return pgaColor.at(index);
}

static QColor getGradientColorfromGal(float value)
{
    //float pg = value / 980 * 100; // convert gal to %g
    float pg = value;

    QList<float> pgaRange;
    //pgaRange << 0 << 0.1 << 0.3 << 0.5 << 2.4 << 6.7 << 13 << 24 << 44 << 83;
    pgaRange << 0 << 0.98 << 2.94 << 4.90 << 23.52 << 65.66 << 127.40 << 235.20 << 431.20 << 813.40;
    QList<QColor> pgaColor;
    pgaColor << QColor("#FFFFFF") << QColor("#A5DDF9") << QColor("#92D050") << QColor("#FFFF00")
             << QColor("#FFC000") << QColor("#FF0000") << QColor("#A32777") << QColor("#632523")
             << QColor("#4C2600") << QColor("#000000");

    int index;

    for(int i=1;i<9;i++)
    {
        if(pg < 0.98)
        {
            index = 0;
            break;
        }
        else if(pg > 813.40)
        {
            index = 9;
            return pgaColor.at(index);
            break;
        }
        else
        {
            if(pg >= pgaRange.at(i) && pg < pgaRange.at(i+1))
            {
                index = i;
                break;
            }
        }
    }

    double k = pgaRange.at(index+1) - pgaRange.at(index);
    double k2 = pgaRange.at(index+1) - pg;
    float key = 100 - (k2 / k * 100); // get percent
    float ratio = key / 100;

    QColor startC = pgaColor.at(index + 1);
    QColor endC = pgaColor.at(index);

    int r = (int)(ratio*startC.red() + (1-ratio)*endC.red());
    int g = (int)(ratio*startC.green() + (1-ratio)*endC.green());
    int b = (int)(ratio*startC.blue() + (1-ratio)*endC.blue());

    return QColor::fromRgb(r, g, b);
}

#endif // COMMON_H
