#include "writelog.h"

WriteLog::WriteLog()
{

}

void WriteLog::write(QString dir, QString str)
{
    QDateTime now = QDateTime::currentDateTimeUtc();
    QDir logDir(dir);
    QFile file(logDir.path() + "/ObservedPGAtoPNG." + now.toString("yyyyMMdd") + ".log");
    if(file.open(QIODevice::WriteOnly | QIODevice::Append))
    {
        QTextStream stream( &file );
        stream << now.toString("[hh:mm:ss] ") << str << "\n";
        file.close();
    }
}
