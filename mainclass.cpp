#include "mainclass.h"

#include "landXY.h"
#include "boundXY.h"
#include "savepng.h"

#include <map>
#include <set>
#include <vector>

#include <fcntl.h>
#include <sstream>
#include <fstream>

using namespace std;
typedef multimap<int, int> mm_type ;

MainClass::MainClass(QString dir)
{
    workDir = dir;
    log = new WriteLog();

    initProj();
    readEventFile();
    ll2xy(pj_longlat, pj_eqc, event.lon, event.lat, &originX, &originY);

    readStaPGAFile();

    if(!initSDL())
    {
        log->write(workDir, "Can't initialize SDL. Exit.");
        exit(1);
    }

    if(!loadBackgroundImage())
    {
        log->write(workDir, "Can't read background image file. Exit.");
        exit(1);
    }

    if(!loadFontFile())
    {
        log->write(workDir, "Can't read font file. Exit");
        exit(1);
    }
    gfxPrimitivesSetFont(myfont, 9, 18) ;

    //applySurface(0, 0, background, screen);

    points = initMSM(true);
    bpoints = initMSM(false);

    QFuture<void> future = QtConcurrent::map(points, getMapZValue);
    future.waitForFinished();

    QFuture<void> future2 = QtConcurrent::map(bpoints, getMapZValue);
    future2.waitForFinished();

    float maxPGA = 0;
    for(int i=0;i<points.size();i++)
    {
        _POINT py = points.at(i);
        if(py.mapZ > maxPGA)
            maxPGA = py.mapZ;
    }
    /*
    qDebug() << "MAX gal:" << maxPGA;
    qDebug() << "Number of stations:" << staList.count();
    qDebug() << "Number of Filtered stations:" << filterOutStaList.count();
    */

    applySurface(0, 0, background, screen);
    SDL_LockSurface(screen);
    Uint32 *pixels = (Uint32 *)screen->pixels;

    for(int i=0;i<points.size();i++)
    {
        _POINT po = points.at(i);

        int pointxy = po.landY * screen->w + po.landX;

        float myZ = po.mapZ;
        QColor col = getGradientColorfromGal(myZ);
        //QColor col = getColorfromGal(myZ);
        pixels[pointxy] = SDL_MapRGB(screen->format, col.red(), col.green(), col.blue());
        //pixels[pointxy] = SDL_MapRGB(screen->format, redColor(myZ), greenColor(myZ),
        //                             blueColor(myZ));
    }

    /*
    for(int i=0;i<bpoints.size();i++)
    {
        _POINT po = bpoints.at(i);

        int pointxy = po.landY * screen->w + po.landX;

        pixels[pointxy] = SDL_MapRGB(screen->format, 80, 80, 80);
    }
    */

    //qDebug() << "NSTA:" << staList.size();

    SDL_UnlockSurface(screen);

    getEdge(0.98);
    getEdge(2.94);
    getEdge(4.90);
    getEdge(23.52);
    getEdge(65.66);
    getEdge(127.40);
    getEdge(235.20);
    getEdge(431.20);
    getEdge(813.40);

    //qDebug() << xyLandList.size();

    for(int i=0;i<xyLandList.size();i++)
    {
        float gal;
        QString intenS;
        int inten;
        if(i==0) { gal = 0.98; intenS = "I"; inten = 1; }
        else if(i==1) { gal = 2.94; intenS = "II"; inten = 2;}
        else if(i==2) { gal = 4.90; intenS = "III"; inten = 3;}
        else if(i==3) { gal = 23.52; intenS = "IV"; inten = 4;}
        else if(i==4) { gal = 65.66; intenS = "V"; inten = 5;}
        else if(i==5) { gal = 127.40; intenS = "VI"; inten = 6;}
        else if(i==6) { gal = 235.20; intenS = "VII"; inten = 7;}
        else if(i==7) { gal = 431.20; intenS = "VIII"; inten = 8;}
        else if(i==8) { gal = 813.40; intenS = "IX"; inten = 9;}
        remove_duplicate_xyList(i, i+1, gal, intenS, inten);
    }

    draw_epicenter();
    //draw_sta();

    shot = SDL_PNGFormatAlpha(screen);

    QString pngfile = workDir + "/ObservedPGA_" + event.evid + ".png";
    const QByteArray stringData = pngfile.toUtf8();
    char text[80];
    text[qMin(79, stringData.size())]='\0';
    std::copy(stringData.constBegin(), stringData.constBegin()+qMin(79, stringData.size()), text);
    SDL_SavePNG(shot, text);
    SDL_FreeSurface(shot);
    SDL_FreeSurface(background);
    SDL_Quit();

    saveStaListToOutput();
    //save_all_points();
    exit(1);
}

void MainClass::draw_epicenter()
{
    Sint16 sx, sy;
    sx = originX - 8;
    sy = originY - 4;
    Sint16 x[10];
    Sint16 y[10];
    x[0]=sx; x[1]=sx+6; x[2]=sx+8; x[3]=sx+10; x[4]=sx+16;
    x[5]=sx+12; x[6]=sx+14; x[7]=sx+8; x[8]=sx+2; x[9]=sx+4;
    y[0]=sy; y[1]=sy; y[2]=sy-6; y[3]=sy; y[4]=sy;
    y[5]=sy+4; y[6]=sy+10; y[7]=sy+6; y[8]=sy+10; y[9]=sy+4;
    filledPolygonRGBA(screen, x, y, 10, 255, 0, 0, 255);
}

QList<_POINT> MainClass::initMSM(bool isLand)
{
    int i, j;
    QList<_POINT> tempPoints;
    double maxdist = 0;
    double Rw = 0, Rp = 0;
    int Np = 9;
    int Nw = 18;

    /*
    if(isLand)
    {
        float myMaxPGA = 0; // get max station

        for(i=0;i<staList.size();i++)
        {
            _STATION sta = staList.at(i);
            if(myMaxPGA < sta.pga)
            {
                maxX = sta.mapX;
                maxY = sta.mapY;
                myMaxPGA = sta.pga;
            }
        }

        qDebug() << myMaxPGA;
    }
    */

    //maxdist = MAX_DISTANCE;
    for(i=0;i<staList.size();i++)
    {
        for(j=0;j<staList.size();j++)
        {
            _STATION sta1 = staList.at(i);
            _STATION sta2 = staList.at(j);
            double d = sqrt(pow((double)(sta1.mapX - sta2.mapX), 2) +
                            pow((double)(sta1.mapY - sta2.mapY), 2));
            if(maxdist < d)
                maxdist = d;
        }
    }

    Rw = maxdist / 2 * sqrt((double)Nw / (double)staList.size());
    Rp = maxdist / 2 * sqrt((double)Np / (double)staList.size());

    double distan, dis, W;

    int cnt = 0;
    if(isLand)
        cnt = LANDXYCNT;
    else
        cnt = BOUNDXYCNT;

    for(i=0;i<cnt;i++)
    {
        _POINT mypoint;
        mypoint.index = i;
        if(isLand)
        {
            mypoint.landX = landXY[i][0];
            mypoint.landY = landXY[i][1];
        }
        else
        {
            mypoint.landX = boundXY[i][0];
            mypoint.landY = boundXY[i][1];
        }

        mypoint.staList = staList;

        for(j=0;j<staList.size();j++)
        {
            _STATION mysta = staList.at(j);
            distan = sqrt(pow((double)(mypoint.landX - mysta.mapX), 2) +
                          pow((double)(mypoint.landY - mysta.mapY), 2) );
            dis = Rp - distan;
            W = 0;

            if(dis >=0)
                W = pow((dis/(Rp*distan)), 2);

            mypoint.mapLUT.append(W);
        }

        double dist = sqrt(pow((double)(mypoint.landX - originX), 2) +
                           pow((double)(mypoint.landY - originY), 2) );

        mypoint.predPGA = getPredictedValue(dist, event.ml, EQ_DEPTH);

        tempPoints.append(mypoint);
    }

    return tempPoints;
}

void MainClass::saveStaListToOutput()
{
    // orderdFilterOutStaList
    for(int i=0;i<staList.size();i++)
    {
        _STATION sta = staList.at(i);
        if(sta.staName.startsWith("VIRTUAL"))
            continue;

        if(orderdFilterOutStaList.isEmpty())
        {
            orderdFilterOutStaList.push_back(sta);
            continue;
        }

        int isInserted = 0;

        for(int j=0;j<orderdFilterOutStaList.size();j++)
        {
            _STATION sta2 = orderdFilterOutStaList.at(j);
            if(sta.pga >= sta2.pga)
            {
                orderdFilterOutStaList.insert(j, sta);
                isInserted = 1;
                break;
            }
        }

        if(isInserted == 0)
            orderdFilterOutStaList.push_back(sta);
    }

    /*
    for(int i=0;i<orderdFilterOutStaList.size();i++)
    {
        _STATION sta = orderdFilterOutStaList.at(i);
        qDebug() << sta.staName << sta.pga;
    }
    */

    QFile file(workDir + "/staPGAList_" + event.evid + ".txt");

    if(file.open(QIODevice::WriteOnly))
    {
        QTextStream stream( &file );
        stream << "Used " << QString::number(staList.count()-1) << " stations.\n";
        file.close();
    }
    if(file.open(QIODevice::WriteOnly | QIODevice::Append))
    {
        QTextStream stream( &file );
        for(int i=0;i<orderdFilterOutStaList.size();i++)
        {
            _STATION sta = orderdFilterOutStaList.at(i);
            stream << sta.staName << " " << QString::number(sta.pga, 'f', 6)
                   << " " << QString::number(sta.pga / 9.8, 'f', 2) << "\n";
        }
        file.close();
    }
}

void MainClass::draw_sta()
{
    for(int i=0;i<staList.size();i++)
    {
        _STATION sta = staList.at(i);
        filledCircleRGBA(screen, sta.mapX, sta.mapY, 1, 0, 0, 255, 255);
/*
        if(sta.pga > 0)
            filledCircleRGBA(screen, sta.mapX, sta.mapY, 1, 0, 0, 255, 255);
        else
            filledCircleRGBA(screen, sta.mapX, sta.mapY, 1, 255, 0, 0, 255);
*/
    }

    for(int i=0;i<filterOutStaList.size();i++)
    {
        _STATION sta = filterOutStaList.at(i);
        filledCircleRGBA(screen, sta.mapX, sta.mapY, 1, 255, 0, 0, 255);
        qDebug() << sta.staName;
    }
}

void MainClass::ll2xy(projPJ src, projPJ target, float lon, float lat, int *rx, int *ry)
{
    double x, y;
    int rtn = 0;
    double mapx, mapy;

    x = lon * DEG_TO_RAD;
    y = lat * DEG_TO_RAD;

    rtn = pj_transform(src, target, 1, 1, &x, &y, nullptr);
    mapx = x - XORIGIN;
    mapx *= SCREEN_WIDTH;
    mapx /= XWIDTH;

    mapy = y - YORIGIN;
    mapy *= SCREEN_HEIGHT;
    mapy /= YHEIGHT;
    mapy = SCREEN_HEIGHT - mapy;

    *rx = (int)mapx;
    *ry = (int)mapy;
}

void MainClass::applySurface(int x, int y, SDL_Surface *source, SDL_Surface *destination)
{
    SDL_Rect offset;
    offset.x = x;
    offset.y = y;
    SDL_BlitSurface(source, nullptr, destination, &offset);
}

void MainClass::initProj()
{
    if(!(pj_longlat = pj_init_plus("+proj=longlat +ellps=WGS84")) )
    {
        qDebug() << "Can't initialize projection.";
        exit(1);
    }
    if(!(pj_eqc = pj_init_plus("+proj=eqc +ellps=WGS84")) )
    {
        qDebug() << "Can't initialize projection.";
        exit(1);
    }
}

bool MainClass::initSDL()
{
    if(SDL_Init(SDL_INIT_VIDEO) == -1)
    {
        return false;
    }

    screen = SDL_CreateRGBSurface( SDL_SWSURFACE, SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, 0, 0, 0, 0);

    if(screen == nullptr)
    {
        return false;
    }

    SDL_WM_SetCaption("Observed PGA. K-SEIS 2020.....KIGAM", NULL);

    c_black = SDL_MapRGB(screen->format,   0,   0,   0);
    c_red   = SDL_MapRGB(screen->format, 255,   0,   0);
    c_white = SDL_MapRGB(screen->format, 255, 255, 255);
    c_gray  = SDL_MapRGB(screen->format, 200, 200, 200);
    c_dgray = SDL_MapRGB(screen->format,  64,  64,  64);
    c_cyan  = SDL_MapRGB(screen->format,  32, 255, 255);

    return true;
}

bool MainClass::loadBackgroundImage()
{
    //background = load_image("/data/4_GIT/newObservedPGAtoPNG/newKorea.png");
    background = load_image("/home/sysop/data/4_GIT/KGShakeMap/V1_0/skorea_2020_07_09.png");

    if(background == nullptr)
    {
        return false;
    }
    return true;
}

SDL_Surface *MainClass::load_image(std::string filename)
{
    SDL_Surface* loadedImage = nullptr;
    loadedImage = IMG_Load(filename.c_str());

    return loadedImage;
}

void MainClass::readStaPGAFile()
{
    double maxPGA = 0;
    QString maxSTA;

    QFile staFile(workDir + "/staPGAList.txt");
    if(!staFile.exists())
    {
        log->write(workDir, "The staList.txt file doesn't exists.");
        exit(1);
    }

    if(staFile.open(QIODevice::ReadOnly))
    {
        QTextStream stream(&staFile);
        QString line, _line;
        int staIndex;

        while(!stream.atEnd())
        {
            line = stream.readLine();
            _line = line.simplified();

            _STATION sta;
            sta.staName = _line.section(" ", 0, 0);

            staIndex = staNameList.indexOf(sta.staName);

            if(staIndex != -1)
                continue;
            else
                staNameList.append(sta.staName);

            sta.lat = _line.section(" ", 1, 1).toDouble();
            sta.lon = _line.section(" ", 2, 2).toDouble();
            ll2xy(pj_longlat, pj_eqc, sta.lon, sta.lat, &sta.mapX, &sta.mapY);
            sta.pga = _line.section(" ", 3, 3).toFloat() * 980;
            //sta.pga = _line.section(" ", 3, 3).toFloat();
            sta.distance = getDistance(sta.lat, sta.lon, event.lat, event.lon);
            sta.predPGA = getPredictedValue(sta.distance, event.ml, EQ_DEPTH);

            //qDebug() << sta.staName << sta.distance << sta.pga << sta.predPGA;

            if(sta.pga >= sta.predPGA - (sta.predPGA * VALID_RANGE) &&
                    sta.pga <= sta.predPGA + (sta.predPGA * VALID_RANGE))
            {
                staList.append(sta);
                if(maxPGA < sta.pga)
                {
                    maxPGA = sta.pga;
                    maxSTA = sta.staName;
                }
            }
            else
            {
                filterOutStaList.append(sta);
            }
        }
        staFile.close();
    }

    //qDebug() << maxPGA << maxSTA;

    // adding one virtual station on origin
    _STATION osta;
    osta.lat = event.lat;
    osta.lon = event.lon;
    osta.mapX = originX;
    osta.mapY = originY;
    osta.pga = maxPGA;
    //osta.pga = getPredictedValue(0, event.ml, EQ_DEPTH);
    osta.staName = "VIRTUAL";
    staList.append(osta);
}

void MainClass::readEventFile()
{
    //<earthquake id="KG21506" netid="KG" network="KIGAM KISS NETWORK" lat="35.3908" lon="128.0092" depth="8" mag="2.8" time="
    //2020-03-20T07:02:46Z" locstring="경남 산청군 동남동쪽 약 13km" event_type="ACTUAL"/

    QFile eventFile(workDir + "/event.xml");
    if(!eventFile.exists())
    {
        log->write(workDir, "The event.xml file doesn't exists.");
        exit(1);
    }

    if(eventFile.open(QIODevice::ReadOnly))
    {
        QTextStream stream(&eventFile);
        QString line, _line;
        line = stream.readLine();
        _line = line.simplified();
        _line.replace(QString("\n"), QString(""));

        QXmlStreamReader xml(_line);
        while(!xml.atEnd() && !xml.hasError())
        {
            xml.readNext();
            if(xml.name().startsWith("earthquake"))
            {
                QXmlStreamAttributes attributes = xml.attributes();
                event.evid = attributes.value("id").toString();
                event.lat = attributes.value("lat").toDouble();
                event.lon = attributes.value("lon").toDouble();
                event.depth = attributes.value("depth").toDouble();
                event.loc = attributes.value("locstring").toString();
                event.origintime = QDateTime::fromString(attributes.value("time").toString(), "yyyy-MM-ddThh:mm:ssZ");
                event.ml = attributes.value("mag").toDouble();
                break;
            }
        }

        eventFile.close();
    }

    log->write(workDir, "EVID:" + event.evid);
    log->write(workDir, "LAT:" + QString::number(event.lat, 'f', 4));
    log->write(workDir, "LON:" + QString::number(event.lon, 'f', 4));
    log->write(workDir, "DEPTH:" + QString::number(event.depth, 'f', 2) + "km");
    log->write(workDir, "LOC:" + event.loc);
    log->write(workDir, "ORIGINTIME:" + event.origintime.toString("yyyy-MM-dd hh:mm:ss"));
    log->write(workDir, "MAG:" + QString::number(event.ml, 'f', 1));
}

void MainClass::getEdge(float gal)
{
    mm_type mm;
    set< pair<int, int> > myequiv;

    set<int> uniq_keys ;

    QList<_POINT> tempPoints;
    for(int i=0;i<points.size();i++)
    {
        _POINT myp = points.at(i);
        if(myp.mapZ >= gal)
        {
            tempPoints.append(myp);
        }
    }
    for(int i=0;i<bpoints.size();i++)
    {
        _POINT myp = bpoints.at(i);
        if(myp.mapZ >= gal)
        {
            tempPoints.append(myp);
        }
    }

    if(tempPoints.isEmpty())
        return;

    // find min, max for each axis
    int xmin = 9999 ;
    int ymin = 9999 ;
    int xmax = 0 ;
    int ymax = 0 ;

    int find_uniq=0 ;
    for (int i=0; i < tempPoints.size(); i ++ )
    {
        _POINT myp = tempPoints.at(i);
        if ( xmin > myp.landX ) xmin = myp.landX ;
        if ( ymin > myp.landY ) ymin = myp.landY ;
        if ( xmax < myp.landX ) xmax = myp.landX ;
        if ( ymax < myp.landY ) ymax = myp.landY ;

        mm.insert(pair<int, int>(myp.landY, myp.landX));
        if( find_uniq != myp.landY )
        {
            uniq_keys.insert(myp.landY) ;
            find_uniq = myp.landY ;
        }
    }
    // make 2d array
    vector< vector<int> > array2D;

    // Set up sizes. (HEIGHT x WIDTH) with margin +- 2
    int width = xmax - xmin + 4 ;
    int height = ymax - ymin + 4 ;

    array2D.resize(height);
    for (int i = 0; i < height; ++i)
      array2D[i].resize(width);

    for (int i = 0; i < height; ++i)
      for ( int j = 0 ; j < width ; ++j)
        array2D[i][j] = 0 ;

    set<int>::reverse_iterator iter;
    mm_type::iterator miter ;
    mm_type::iterator lower_mit ;
    mm_type::iterator upper_mit ;

    int label = 0 ;
    int col, row ;

    // Connected Component Labelling algorithm
    for (iter = uniq_keys.rbegin(); iter != uniq_keys.rend(); ++iter)
    {
        lower_mit = mm.lower_bound(*iter) ;
        upper_mit = mm.upper_bound(*iter) ;

        // width
        for ( miter=lower_mit; miter!=upper_mit; ++miter)
        {
            col = (*miter).second - xmin + 2 ;
            row = (*miter).first - ymin + 2 ;

            if ( array2D[row+1][col] == 0 && array2D[row][col-1] == 0 )
            {
                array2D[row][col] = ++ label ;
            }
            else if ( array2D[row+1][col] != 0 && array2D[row][col-1] == 0 )
            {
                array2D[row][col] = array2D[row+1][col] ;
            }
            else if ( array2D[row][col-1] != 0 && array2D[row+1][col] == 0 )
            {
                array2D[row][col] = array2D[row][col-1] ;
            }
            else if ( array2D[row][col-1] != 0 && array2D[row+1][col] != 0 )
            {
                array2D[row][col] = array2D[row+1][col] ;
                if ( array2D[row][col-1] != array2D[row+1][col] )
                {
                    myequiv.insert(pair<int, int>(array2D[row+1][col], array2D[row][col-1]));
                }
            }
        }
    }

    for (set< pair<int, int> >::reverse_iterator mi = myequiv.rbegin(); mi != myequiv.rend(); ++mi)
    {
        for (int i = 0; i < height; ++i)
        {
            for ( int j = 0; j < width; ++j)
            {
                if ( array2D[i][j] == (*mi).second )
                    array2D[i][j] = (*mi).first ;
            }
        }
    }

    // now edge detect
    for (iter = uniq_keys.rbegin(); iter != uniq_keys.rend(); ++iter)
    {
        lower_mit = mm.lower_bound(*iter) ;
        upper_mit = mm.upper_bound(*iter) ;

        // width
        for ( miter=lower_mit; miter!=upper_mit; ++miter)
        {
            col = (*miter).second - xmin + 2 ;
            row = (*miter).first - ymin + 2 ;

            if ( array2D[row+1][col] == 0 || array2D[row-1][col] == 0 ||
                 array2D[row][col-1] == 0 || array2D[row][col+1] == 0 )
                array2D[row][col] = 1 ;
            else
                array2D[row][col] = 2 ;
        }
    }

    // now set to 0 except edge
    for (int i=0; i<height;++i)
    {
        for (int j=0;j<width;++j)
        {
            if ( array2D[i][j] == 2 )
            {
                array2D[i][j] = 0 ;
            }
        }
    }

    XYLAND xyland;

    for (int i = 0; i < height; ++i)
    {
        for ( int j = 0; j < width; ++j)
        {
            if(array2D[i][j])
            {
                filledCircleRGBA(screen, j+xmin-2, i+ymin-2, 0.2, 0, 0, 0, 255);

                xyland.xList.append(j+xmin-2);
                xyland.yList.append(i+ymin-2);
            }
        }
    }

    xyLandList.append(xyland);
}

void MainClass::remove_duplicate_xyList(int first, int second,
                                        float gal, QString intenS, int inten)
{
    XYLAND fXY = xyLandList.at(first);
    int boxx, boxy;
    if(second >= xyLandList.size())
    {
        boxx = fXY.xList.at(fXY.xList.size()/2);
        boxy = fXY.yList.at(fXY.yList.size()/2);
    }
    else
    {
        XYLAND sXY = xyLandList.at(second);
        if(fXY.xList.first() != sXY.xList.first())
        {
            if(fXY.yList.first() <= 15)
            {
                boxx = fXY.xList.at(fXY.xList.size()/2);
                boxy = fXY.yList.at(fXY.yList.size()/2);
            }
            else
            {
                boxx = fXY.xList.first();
                boxy = fXY.yList.first();
            }
        }
        else
        {
            for(int i=0;i<fXY.xList.size();i++)
            {
                int fx = fXY.xList.at(i);
                int fy = fXY.yList.at(i);
                for(int j=0;j<sXY.xList.size();j++)
                {
                    int sx = sXY.xList.at(j);
                    int sy = sXY.yList.at(j);


                    if(fx == sx && fy == sy)
                    {
                        fXY.xList.replace(i, 0);
                        fXY.yList.replace(i, 0);
                        break;
                    }
                }
            }

            fXY.xList.removeAll(0);
            fXY.yList.removeAll(0);

            boxx = fXY.xList.at(fXY.xList.size()/2);
            boxy = fXY.yList.at(fXY.yList.size()/2);
        }
    }

    boxx = boxx - 12;
    boxy = boxy - 8;
    SDL_Rect target_rect;
    target_rect.x = boxx;
    target_rect.y = boxy;

    int boxw;

    /*
    if(inten == 1) boxw = 49;
    else if(inten == 2) boxw = 56;
    else if(inten == 3) boxw = 65;
    else if(inten == 4) boxw = 57;
    else if(inten == 5) boxw = 50;
    else if(inten == 6) boxw = 50;
    else if(inten == 7) boxw = 56;
    else if(inten == 8) boxw = 64;
    else boxw = 53;
    */

    boxw = 31;

    target_rect.w = boxw;
    target_rect.h = 15;
    SDL_FillRect(screen, &target_rect, c_white);
    rectangleRGBA( screen, target_rect.x, target_rect.y,
                   target_rect.x+boxw, target_rect.y+15 , 0, 0, 0, 255);

    char tstr[80] ={0x00,};

    gal = gal / 980 * 100; // convert gal to %g

    QString mystr;

    /*
    if(inten >= 6)
        mystr = QString::number(gal) + "/" + intenS;
    else
        mystr = QString::number(gal, 'f', 1) + "/" + intenS;
        */

    if(inten >= 6)
        mystr = QString::number(gal);
    else
        mystr = QString::number(gal, 'f', 1);

    memcpy(tstr, mystr.toStdString().c_str(), mystr.size());

    stringRGBA(screen, boxx + 2, boxy, tstr, 0, 0, 0, 255);

}

void MainClass::draw_legend()
{
    char tstr[80] ={0x00,} ;

    SDL_Rect r1; r1.x = 297; r1.y=708; r1.w=63; r1.h=26;
    SDL_FillRect(screen, &r1, c_gray);
    rectangleRGBA( screen, r1.x, r1.y,
                   r1.x+63, r1.y+26 , 0, 0, 0, 255);

    sprintf(tstr, "PGA(%%g)") ;
    stringRGBA(screen, r1.x+2, r1.y+4, tstr, 0, 0, 0, 255);

    SDL_Rect r2; r2.x = 297; r2.y=708+26; r2.w=63; r2.h=26;
    SDL_FillRect(screen, &r2, c_gray);
    rectangleRGBA( screen, r2.x, r2.y,
                   r2.x+63, r2.y+26 , 0, 0, 0, 255);

    sprintf(tstr, "MMI") ;
    stringRGBA(screen, r2.x+20, r2.y+4, tstr, 0, 0, 0, 255);

    int bwidth = 45;
    int bheight = 26;
    int xstart = 360;
    int ystart = 708;

    float mygal;
    QColor col;
    Uint32 myc;
    QString str;

    QList<float> pgaRange;
    pgaRange << 0 << 0.98 << 2.94 << 4.90 << 23.52 << 65.66 << 127.40 << 235.20 << 431.20 << 813.40;

    for(int i=0;i<10;i++)
    {
        SDL_Rect t1;
        t1.x = xstart + (i*bwidth);
        t1.y = ystart;
        t1.w = bwidth ;
        t1.h = bheight ;

        for(int j=0;j<9;j++)
        {
            if(i!=9)
            {
                float mygal = pgaRange[i] + ((pgaRange[i+1] - pgaRange[i])/9*j);
                QColor col = getGradientColorfromGal(mygal);
                myc  = SDL_MapRGB(screen->format, col.red(), col.green(), col.blue());
            }
            else if(i==9)
                myc = SDL_MapRGB(screen->format, 0, 0, 0);

            SDL_Rect t2;
            t2.x = t1.x + (j*5);
            t2.y = t1.y;
            t2.w = 5;
            t2.h = bheight;
            SDL_FillRect(screen, &t2, myc);
        }
        rectangleRGBA( screen, t1.x, t1.y,
                       t1.x+bwidth, t1.y+bheight , 0, 0, 0, 255);

        setlocale(LC_ALL, "");
        //int r = wctomb(tstr, 0x2264);
        //tstr[r] = '\0';

        if(i==0) sprintf(tstr, "<0.1");
        else if(i==1) sprintf(tstr, "%lc 0.3", L'≤');
        else if(i==2) sprintf(tstr, "%lc 0.5", L'≤');
        else if(i==3) sprintf(tstr, "%lc 2.4", L'≤');
        else if(i==4) sprintf(tstr, "\xF3 ≤6.7");
        else if(i==5) sprintf(tstr, "\342\211\244 13");
        else if(i==6) sprintf(tstr, "\u2264");
        else if(i==7) sprintf(tstr, "≤44");
        else if(i==8) sprintf(tstr, "≤83");
        else if(i==9) sprintf(tstr, ">83");

       // wprintf(L"Here are ≤");

        //QChar ttt[1] = "≤";
        //qDebug() << ttt.toLatin1();

        stringRGBA(screen, t1.x+3, t1.y+4, tstr, 0, 0, 0, 255);

        t1.y = ystart + bheight ;
        for(int j=0;j<9;j++)
        {
            if(i!=9)
            {
                float mygal = pgaRange[i] + ((pgaRange[i+1] - pgaRange[i])/9*j);
                QColor col = getGradientColorfromGal(mygal);
                myc  = SDL_MapRGB(screen->format, col.red(), col.green(), col.blue());
            }
            else if(i==9)
                myc = SDL_MapRGB(screen->format, 0, 0, 0);
            SDL_Rect t2;
            t2.x = t1.x + (j*5);
            t2.y = t1.y;
            t2.w = 5;
            t2.h = bheight;
            SDL_FillRect(screen, &t2, myc);
        }
        rectangleRGBA( screen, t1.x, t1.y,
                       t1.x+bwidth, t1.y+bheight , 0, 0, 0, 255);
    }
}

bool MainClass::loadFontFile()
{
    char filename[128] = "/data/4_GIT/ObservedPGAtoPNG/9x18B.fnt";
    //char filename[128] =
    //        "/data/4_GIT/ObservedPGAtoPNG/Cantarell-Regular.otf";
    int filelen = 9216 ;
    FILE *file;
    int bytesRead;

    /* check file exist */
    if(std::ifstream(filename))
    {
        /* Allocate memory for font data */
        myfont=(char *)malloc(filelen) ;
        if (myfont)
        {
            file = fopen(filename,"r");
            bytesRead = fread(myfont,filelen,1,file);
            fclose(file);
            return true ;
        }
    }
    else
        return false ;
}

void MainClass::xy2ll(projPJ src, projPJ target,
                      int x, int y, float gal)
{
    double xx, yy;
    int rtn = 0;
    xx = x * XWIDTH;
    xx = xx / SCREEN_WIDTH;
    xx = xx + XORIGIN;

    //yy = SCREEN_HEIGHT + y;
    yy = SCREEN_HEIGHT - y;
    yy = yy * YHEIGHT;
    yy = yy / SCREEN_HEIGHT;
    yy = yy + YORIGIN;

    rtn = pj_transform(src, target, 1, 1, &xx, &yy, nullptr);

    xx = xx / DEG_TO_RAD;
    yy = yy / DEG_TO_RAD;

    QFile file(workDir + "/GALforAllPoints_" + event.evid + ".txt");
    if(file.open(QIODevice::WriteOnly | QIODevice::Append))
    {
        QTextStream stream( &file );
        stream << QString::number(xx, 'f', 4) << " "
               << QString::number(yy, 'f', 4) << " "
               << QString::number(gal, 'f', 4) << "\n";
        file.close();
    }
}

void MainClass::save_all_points()
{

    for(int i=0;i<points.size();i++)
    {
        _POINT po = points.at(i);
        xy2ll(pj_eqc, pj_longlat, po.landX, po.landY, po.mapZ);
    }

    for(int i=0;i<bpoints.size();i++)
    {
        _POINT po = bpoints.at(i);
        xy2ll(pj_eqc, pj_longlat, po.landX, po.landY, po.mapZ);
    }
}
