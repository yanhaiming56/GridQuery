#ifndef COMALGORITHM_H
#define COMALGORITHM_H

#include "ClassType.h"
#include "../ine_include/INE.h"

using namespace std;

class ComAlgorithm {
public:
    ComAlgorithm();

    virtual ~ComAlgorithm();

public:
    int getGridNum() {
        return m_gridNum;
    }

    void setGridNum(int gridNum) {
        m_gridNum = gridNum;
    }

    double getGridUnit() {
        return m_gridUnit;
    }

    void setGridUnit(double gridUnit) {
        m_gridUnit = gridUnit;
    }

    //获取外部数据
public:
    void getVertexs(const string &strFileName, Vertexs *&vertexs);

    void getDataPoints(const string &strFileName, DataPoints *&dataPoints);

    void getRoads(const string &strFileName, Vertexs *const &vertexs, Roads *&roads);

    void getRoadOfPoints(const string &strFileName, DataPoints *&dataPoints, Roads *&roads, Vertexs *const &vertexs);


    void getNearestPoint(const string &strFileName, Roads *const &roads, Vertexs *const &vertexs,
                         DataPoints *const &dataPoints, INE &ine);

    void statisticSubRoadofGtids(Roads *const &roads, Vertexs *const &vertexs,
                                 DataPoints *const &dataPoints);

    string &trimString(string &s);

    void createQueryData(Roads *const &roads);

    DataPoint getNearestDataOfVertexs(Vertexs *const &vertexs, const int &vid, Roads *const &roads, const int &rid,
                                      DataPoints *const &dataPoints, INE &ine);

    void NaiveNN(const string &strFileName, Roads *const &roads, DataPoints *const &datapoints);

    int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

private:
    Point getMidstOfTwoPoint(const Point &startPoint, const Point &endPoint);

    Point calPointFromDistance(const double &distance, const Road &road);

    void subRoadPassGrids(const SubRoad &subRoad);

    double getSlope(const Point &startPoint, const Point &endPoint);

    double getDistanceOfTwoPoint(const Point &startPoint, const Point &endPoint);


protected:

private:
    size_t m_gridNum;
    double m_gridUnit;
    Point m_MinPoint;
    Point m_MaxPoint;
    Grids *m_Grids;

};

#endif // COMALGORITHM_H
