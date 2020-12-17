#ifndef MAINCLASS_H
#define MAINCLASS_H

#include <QObject>
#include <QXmlStreamReader>
#include <QXmlStreamAttributes>
#include <QMap>
#include <QtConcurrent>

#include "common.h"
#include "writelog.h"

#include "SDL/SDL.h"
#include "SDL/SDL_image.h"
#include "SDL/SDL_gfxPrimitives.h"

#include "wchar.h"

class MainClass
{
public:
    MainClass(QString dir = nullptr);

private:
    WriteLog *log;

    void initProj();
    projPJ pj_eqc;
    projPJ pj_longlat;
    void ll2xy(projPJ src, projPJ target, float lon, float lat, int *rx, int *ry);
    void xy2ll(projPJ src, projPJ target,
               int x, int y, float gal);

    bool initSDL();
    SDL_Surface *background = nullptr;
    SDL_Surface *screen = nullptr;
    SDL_Surface *shot = nullptr;
    Uint32 c_black;
    Uint32 c_red;
    Uint32 c_white;
    Uint32 c_gray;
    Uint32 c_dgray;
    Uint32 c_cyan;
    bool loadBackgroundImage();
    void applySurface(int x, int y, SDL_Surface* source, SDL_Surface* destination);
    SDL_Surface *load_image( std::string filename );

    QString workDir;
    _EVENT event;
    QList<_STATION> staList;
    QList<_STATION> filterOutStaList;
    QList<_STATION> orderdFilterOutStaList;
    QList<QString> staNameList;

    QList<_POINT> points;
    QList<_POINT> bpoints;

    void readEventFile();
    void readStaPGAFile();

    QList<_POINT> initMSM(bool type);
    void getEdge(float);

    int originX, originY;
    int maxX, maxY;

    QList<XYLAND> xyLandList;

    char *myfont ;

    void saveStaListToOutput();
    void draw_sta();
    void draw_legend();
    bool loadFontFile();
    void draw_epicenter();
    void remove_duplicate_xyList(int, int, float, QString, int);
    void save_all_points();
};

#endif // MAINCLASS_H
