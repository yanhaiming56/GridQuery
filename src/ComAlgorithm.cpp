#include "ComAlgorithm.h"
#include "INEClassType.h"

using namespace std;

ComAlgorithm::ComAlgorithm()
        : m_gridNum(0), m_gridUnit(0), m_Grids(nullptr) {
    //ctor
}

ComAlgorithm::~ComAlgorithm() {
    //dtor
    if (m_Grids != nullptr) {
        for (size_t i = 0; i < m_gridNum; i++)
            delete[] m_Grids[i].m_Grids;

        delete[] m_Grids;
    }
}

/** @brief getVertexs

    @todo: document this function
*/
void ComAlgorithm::getVertexs(const string &strFileName, Vertexs *&vertexs) {
    if (!strFileName.compare("")) {
        cout << "顶点文件错误！" << endl;
        return;
    }

    if (vertexs == nullptr)
        vertexs = new Vertexs[MAXLEN];

    ifstream ifile(strFileName);
    string strBuffer;

    double startMin = MAXNUM, endMin = MAXNUM;
    double startMax = MINNUM, endMax = MINNUM;
    int lineNum = -1;

    while (getline(ifile, strBuffer)) {
        strBuffer = trimString(strBuffer);

        if (strBuffer.empty())
            continue;

        vector<string> buffers;
        boost::split(buffers, strBuffer, boost::is_any_of(" "));
        Point point(stod(buffers[1]), stod(buffers[2]));

        if (startMin > point.m_XPos)
            startMin = point.m_XPos;

        if (endMin > point.m_YPos)
            endMin = point.m_YPos;

        if (startMax < point.m_XPos)
            startMax = point.m_XPos;

        if (endMax < point.m_YPos)
            endMax = point.m_YPos;

        lineNum = stoi(buffers[0]);

        if (lineNum % MAXLEN == 0) //防止编号从1开始没有开辟空间
            vertexs[lineNum / MAXLEN].m_Vertexs = new Vertex[MAXLEN];

        Vertex v = vertexs[lineNum / MAXLEN].m_Vertexs[lineNum % MAXLEN] = Vertex(lineNum, point);
    }

    m_MinPoint = Point(startMin, endMin);
    m_MaxPoint = Point(startMax, endMax);
    m_gridNum = int(double(startMax - startMin) / double(m_gridUnit)) + 10;
    ifile.close();

}

/** @brief getDataPoints

    @todo: document this function
*/
void ComAlgorithm::getDataPoints(const string &strFileName, DataPoints *&dataPoints) {
    if (!strFileName.compare("")) {
        cout << "数据点文件错误！" << endl;
        return;
    }

    if (dataPoints == nullptr)
        dataPoints = new DataPoints[MAXLEN];

    ifstream ifile(strFileName);
    string strBuffer;

    int lineNum = -1;

    while (getline(ifile, strBuffer)) {
        strBuffer = trimString(strBuffer);

        if (strBuffer.empty())
            continue;

        vector<string> buffers;
        boost::split(buffers, strBuffer, boost::is_any_of(" "));
        Point point(stod(buffers[1]), stod(buffers[2]));

        lineNum = stoi(buffers[0]);

        if (lineNum % MAXLEN == 0)
            dataPoints[lineNum / MAXLEN].m_DataPoints = new DataPoint[MAXLEN];

        dataPoints[lineNum / MAXLEN].m_DataPoints[lineNum % MAXLEN] = DataPoint(lineNum, 0, point);
    }
}

/** @brief getRoads

    @todo: document this function
*/
void ComAlgorithm::getRoads(const string &strFileName, Vertexs *const &vertexs, Roads *&roads) {
    if (!strFileName.compare("")) {
        cout << "道路文件错误！" << endl;
        return;
    }

    if (vertexs == nullptr) {
        cout << "顶点数据错误！" << endl;
        return;
    }

    if (roads == nullptr)
        roads = new Roads[MAXLEN];

    ifstream ifile(strFileName);
    string strBuffer;

    while (getline(ifile, strBuffer)) {
        vector<string> buffers;
        boost::split(buffers, strBuffer, boost::is_any_of(" "));

        int lineNum = stoi(buffers[0]);

        if (lineNum % MAXLEN == 0)
            roads[lineNum / MAXLEN].m_Roads = new Road[MAXLEN];

        Vertex *vStart = &(vertexs[stoi(buffers[1]) / MAXLEN].m_Vertexs[stoi(
                buffers[1]) % MAXLEN]);
        vStart->m_Rids.push_back(lineNum);
        Vertex *vEnd = &(vertexs[stoi(buffers[2]) / MAXLEN].m_Vertexs[stoi(buffers[2]) % MAXLEN]);
        vEnd->m_Rids.push_back(lineNum);

        roads[lineNum / MAXLEN].m_Roads[lineNum % MAXLEN] = Road(lineNum, stod(buffers[3]), *vStart,
                                                                 *vEnd);
    }
}

/** @brief getRoadOfPoints

    @todo: document this function
*/
void ComAlgorithm::getRoadOfPoints(const string &strFileName, DataPoints *&dataPoints,
                                   Roads *&roads, Vertexs *const &vertexs) {
    if (!strFileName.compare("")) {
        cout << "道路数据点关系文件错误！" << endl;
        return;
    }

    if (dataPoints == nullptr)
        dataPoints = new DataPoints[MAXLEN];

    if (roads == nullptr)
        roads = new Roads[MAXLEN];

    ifstream ifile(strFileName);
    string strBuffer;

    int pointNum = 0;
    int lineNum = 0;
    bool isRoadFlag = true;
    Road curRoad;

    while (getline(ifile, strBuffer)) {
        strBuffer = trimString(strBuffer);

        if (strBuffer.empty())
            continue;

        vector<string> buffers;
        boost::split(buffers, strBuffer, boost::is_any_of(" "));

        size_t splitCount = buffers.size();
        int numPoints = stoi(buffers[splitCount - 1]);

        if (lineNum % MAXLEN == 0)
            roads[lineNum / MAXLEN].m_Roads = new Road[MAXLEN];

        if (isRoadFlag) {
            Vertex *vStart = &(vertexs[stoi(buffers[0]) / MAXLEN].m_Vertexs[stoi(
                    buffers[0]) % MAXLEN]);
            vStart->m_Rids.push_back(lineNum);
            Vertex *vEnd = &(vertexs[stoi(buffers[1]) / MAXLEN].m_Vertexs[stoi(buffers[1]) % MAXLEN]);
            vEnd->m_Rids.push_back(lineNum);

            curRoad = roads[lineNum / MAXLEN].m_Roads[lineNum % MAXLEN] = Road(lineNum, stod(buffers[2]), *vStart,
                                                                               *vEnd);
            lineNum++;

            if (numPoints != 0)
                isRoadFlag = false;
        } else {
            for (size_t i = 1; i < splitCount; i += 2) {
                double distance = stod(buffers[i]);
                Point dataPoint = calPointFromDistance(distance, curRoad);

                if (pointNum % MAXLEN == 0)
                    dataPoints[pointNum / MAXLEN].m_DataPoints = new DataPoint[MAXLEN];

                DataPoint dt = dataPoints[pointNum / MAXLEN].m_DataPoints[pointNum % MAXLEN] = DataPoint(pointNum,
                                                                                                         distance,
                                                                                                         dataPoint);
                roads[(lineNum - 1) / MAXLEN].m_Roads[(lineNum - 1) % MAXLEN].m_Datas.insert(dt);
                pointNum++;
            }

            isRoadFlag = true;
        }


    }

    ifile.close();
}

Point ComAlgorithm::calPointFromDistance(const double &distance, const Road &road) {
    Point beginPoint = road.m_Begin.m_Point;
    Point endPoint = road.m_End.m_Point;
    /*double slop = road.m_Slope;

    if (fabs(slop) != MAXNUM) {
        double x1_x0 = pow(distance, 2) / (slop * slop + 1);
        double fabs_x1_x0 = sqrt(x1_x0);
        double x1 = beginPoint.m_XPos + fabs_x1_x0;
        double y1 = slop * (x1 - beginPoint.m_XPos) + beginPoint.m_YPos;
        double beginDist = sqrt(pow((x1 - beginPoint.m_XPos), 2) + pow((y1 - beginPoint.m_YPos), 2));
        double endDist = sqrt(pow((x1 - endPoint.m_XPos), 2) + pow((y1 - endPoint.m_YPos), 2));

        if (fabs(beginDist + endDist - road.m_Distance) > ZERO) {
            x1 = beginPoint.m_XPos - fabs_x1_x0;
            y1 = slop * (x1 - beginPoint.m_XPos) + beginPoint.m_YPos;
        }

        return Point(x1, y1);
    } else {
        double x1 = beginPoint.m_XPos;
        double y1 = beginPoint.m_YPos + distance;
        return Point(x1, y1);
    }*/
    double x = 0, y = 0;

    double dx = endPoint.m_XPos - beginPoint.m_XPos;
    double dy = endPoint.m_YPos - beginPoint.m_YPos;
    double b = dx * beginPoint.m_YPos - dy * beginPoint.m_XPos;

    if (fabs(dy) < ZERO) {
        beginPoint.m_XPos > endPoint.m_XPos ? x = beginPoint.m_XPos - distance : x = beginPoint.m_XPos + distance;
        y = beginPoint.m_YPos;
        return Point(x, y);
    }

    if (fabs(dx) < ZERO) {
        beginPoint.m_YPos > endPoint.m_YPos ? y = beginPoint.m_YPos - distance : y = beginPoint.m_YPos + distance;
        x = beginPoint.m_XPos;
        return Point(x, y);
    }
    double slop = dy / dx;
    x = beginPoint.m_XPos + distance / sqrt(1 + slop * slop);
    y = (dy * x + b) / dx;
    double sumLen = getDistanceOfTwoPoint(Point(x, y), beginPoint) + getDistanceOfTwoPoint(Point(x, y), endPoint);
    if (fabs(sumLen - road.m_Distance) < ZERO)
        return Point(x, y);
    else {
        x = beginPoint.m_XPos - distance / sqrt(1 + slop * slop);
        y = (dy * x + b) / dx;
        return Point(x, y);
    }

}

/** @brief getNearestDataOfVertexs

    @todo: document this function
*/
DataPoint
ComAlgorithm::getNearestDataOfVertexs(Vertexs *const &vertexs, const int &vid, Roads *const &roads, const int &rid,
                                      DataPoints *const &dataPoints, INE &ine) {

}


/** @brief statisticSubRoadofGtids

    @todo: document this function
*/
void ComAlgorithm::statisticSubRoadofGtids(Roads *const &roads, Vertexs *const &vertexs,
                                           DataPoints *const &dataPoints) {
    if (m_Grids == nullptr) {
        m_Grids = new Grids[m_gridNum];

        for (size_t i = 0; i < m_gridNum; i++)
            m_Grids[i].m_Grids = new Grid[m_gridNum];
    }

    for (int rowRoad = 0; rowRoad < MAXLEN; rowRoad++) {
        if (roads[rowRoad].m_Roads != nullptr) {
            for (int colRoad = 0; colRoad < MAXLEN; colRoad++) {
                Road road = roads[rowRoad].m_Roads[colRoad];

                if (road.m_Rid == -1)
                    break;

                Vertex v_start = road.m_Begin;
                Vertex v_end = road.m_End;

                set<DataPoint> datasOfRoad = road.m_Datas;

                if (datasOfRoad.size() == 0) {
                    Point midstPoint = getMidstOfTwoPoint(v_start.m_Point, v_end.m_Point);
                    SubRoad startSubRoad(road.m_Rid, v_start.m_Vid, v_start.m_Point, midstPoint, road.m_Slope, true);
                    subRoadPassGrids(startSubRoad);
                    SubRoad endSubRoad(road.m_Rid, v_end.m_Vid, midstPoint, v_end.m_Point, road.m_Slope, true);
                    subRoadPassGrids(endSubRoad);
                }

                if (datasOfRoad.size() == 1) {
                    DataPoint dp = *(datasOfRoad.begin());
                    Point startMidstPoint = getMidstOfTwoPoint(v_start.m_Point, dp.m_Point);
                    SubRoad startSubRoad(road.m_Rid, v_start.m_Vid, v_start.m_Point, startMidstPoint, road.m_Slope,
                                         true);
                    subRoadPassGrids(startSubRoad);

                    Point endMidstPoint = getMidstOfTwoPoint(dp.m_Point, v_end.m_Point);
                    SubRoad midstSubRoad(road.m_Rid, dp.m_Oid, startMidstPoint, endMidstPoint, road.m_Slope, false);
                    subRoadPassGrids(midstSubRoad);

                    SubRoad endSubRoad(road.m_Rid, v_end.m_Vid, endMidstPoint, v_end.m_Point, road.m_Slope, true);
                    subRoadPassGrids(endSubRoad);
                }

                if (datasOfRoad.size() > 1) {
                    DataPoint startData = *(datasOfRoad.begin());
                    Point preMidstPoint = getMidstOfTwoPoint(v_start.m_Point, startData.m_Point);
                    SubRoad startSubRoad(road.m_Rid, v_start.m_Vid, v_start.m_Point, preMidstPoint, road.m_Slope,
                                         true);
                    subRoadPassGrids(startSubRoad);

                    for (set<DataPoint>::iterator it = datasOfRoad.begin(); it != datasOfRoad.end(); it++) {
                        DataPoint curDataPoint = *it;
                        set<DataPoint>::iterator nextIt = ++it;
                        --it;

                        if (nextIt != datasOfRoad.end()) {
                            DataPoint nextDataPoint = *nextIt;
                            Point midstPoint = getMidstOfTwoPoint(curDataPoint.m_Point, nextDataPoint.m_Point);
                            SubRoad curSubRoad(road.m_Rid, curDataPoint.m_Oid, preMidstPoint, midstPoint, road.m_Slope,
                                               false);
                            subRoadPassGrids(curSubRoad);
                            preMidstPoint = midstPoint;

                        }

                        if (nextIt == datasOfRoad.end()) {
                            Point midstPoint = getMidstOfTwoPoint(curDataPoint.m_Point, v_end.m_Point);
                            SubRoad curSubRoad(road.m_Rid, curDataPoint.m_Oid, preMidstPoint, midstPoint, road.m_Slope,
                                               false);
                            subRoadPassGrids(curSubRoad);
                            preMidstPoint = midstPoint;
                        }
                    }

                    SubRoad endSubRoad(road.m_Rid, v_end.m_Vid, preMidstPoint, v_end.m_Point, road.m_Slope, true);
                    subRoadPassGrids(endSubRoad);

                }
            }
        }
    }
}

/** @brief getMidstOfTwoPoint

    @todo: document this function
*/
Point ComAlgorithm::getMidstOfTwoPoint(const Point &startPoint, const Point &endPoint) {
    double xPos = (startPoint.m_XPos + endPoint.m_XPos) / 2;
    double yPos = (startPoint.m_YPos + endPoint.m_YPos) / 2;
    return Point(xPos, yPos);
}

/** @brief subRoadPassGrids

    @todo: document this function
*/
void ComAlgorithm::subRoadPassGrids(const SubRoad &subRoad) {
    Point startPoint = subRoad.m_Start, endPoint = subRoad.m_End;

    if (startPoint.m_XPos > endPoint.m_XPos || (fabs(startPoint.m_XPos - endPoint.m_XPos) < ZERO
                                                && startPoint.m_YPos > endPoint.m_YPos)) {
        Point tempPoint = startPoint;
        startPoint = endPoint;
        endPoint = tempPoint;
    }

    double slope = getSlope(startPoint, endPoint);

    /*
        int x = int ( ( startPoint.m_XPos - m_MinPoint.m_XPos ) / m_gridUnit ) + 1;
        int y = int ( ( startPoint.m_YPos - m_MinPoint.m_YPos ) / m_gridUnit ) + 1;
    */

    if (fabs(slope) >= MAXNUM) {
        int x = int((startPoint.m_XPos - m_MinPoint.m_XPos) / m_gridUnit) + 1;
        int y = int((startPoint.m_YPos - m_MinPoint.m_YPos) / m_gridUnit) + 1;
        Point upper(x * m_gridUnit + m_MinPoint.m_XPos, y * m_gridUnit + m_MinPoint.m_YPos);

        while (upper.m_YPos < endPoint.m_YPos) {
            m_Grids[x].m_Grids[y].insertSubRoad(x, y, subRoad);
            upper = Point(upper.m_XPos, upper.m_YPos + m_gridUnit);
            y++;
        }

        m_Grids[x].m_Grids[y].insertSubRoad(x, y, subRoad);
    } else if (fabs(slope) <= ZERO || fabs(slope) == 0) {
        int x = int((startPoint.m_XPos - m_MinPoint.m_XPos) / m_gridUnit) + 1;
        int y = int((startPoint.m_YPos - m_MinPoint.m_YPos) / m_gridUnit) + 1;

        Point upper(x * m_gridUnit + m_MinPoint.m_XPos, y * m_gridUnit + m_MinPoint.m_YPos);

        while (upper.m_XPos < endPoint.m_XPos) {
            m_Grids[x].m_Grids[y].insertSubRoad(x, y, subRoad);
            upper = Point(upper.m_XPos + m_gridUnit, upper.m_YPos);
            x++;
        }

        m_Grids[x].m_Grids[y].insertSubRoad(x, y, subRoad);
    } else if (slope > 0) {
        int x = int((startPoint.m_XPos - m_MinPoint.m_XPos) / m_gridUnit) + 1;
        int y = int((startPoint.m_YPos - m_MinPoint.m_YPos) / m_gridUnit) + 1;
        Point upper(x * m_gridUnit + m_MinPoint.m_XPos, y * m_gridUnit + m_MinPoint.m_YPos);
        Point preInterPoint = startPoint;

        while (preInterPoint.m_XPos < endPoint.m_XPos) {
            double tempSlope = getSlope(preInterPoint, upper);

            if (slope > tempSlope) {
                double tempX = (upper.m_YPos - preInterPoint.m_YPos) / slope + preInterPoint.m_XPos;
                Point interPoint(tempX, upper.m_YPos);
                double midstX = (preInterPoint.m_XPos + interPoint.m_XPos) / 2;
                double midstY = (preInterPoint.m_YPos + interPoint.m_YPos) / 2;

                int gridX = int((midstX - m_MinPoint.m_XPos) / m_gridUnit) + 1;
                int gridY = int((midstY - m_MinPoint.m_YPos) / m_gridUnit) + 1;

                m_Grids[gridX].m_Grids[gridY].insertSubRoad(gridX, gridY, subRoad);

                double upperY = slope * (upper.m_XPos - startPoint.m_XPos) + startPoint.m_YPos;

                if (fabs(upperY - upper.m_YPos) < ZERO)
                    upper = Point(upper.m_XPos + m_gridUnit, upper.m_YPos + m_gridUnit);
                else
                    upper = Point(upper.m_XPos, upper.m_YPos + m_gridUnit);

                preInterPoint = interPoint;
            } else {
                double tempY = slope * (upper.m_XPos - preInterPoint.m_XPos) + preInterPoint.m_YPos;
                Point interPoint(upper.m_XPos, tempY);
                double midstX = (preInterPoint.m_XPos + interPoint.m_XPos) / 2;
                double midstY = (preInterPoint.m_YPos + interPoint.m_YPos) / 2;

                int gridX = int((midstX - m_MinPoint.m_XPos) / m_gridUnit) + 1;
                int gridY = int((midstY - m_MinPoint.m_YPos) / m_gridUnit) + 1;

                m_Grids[gridX].m_Grids[gridY].insertSubRoad(gridX, gridY, subRoad);

                double upperY = slope * (upper.m_XPos - startPoint.m_XPos) + startPoint.m_YPos;

                if (fabs(upperY - upper.m_YPos) < ZERO)
                    upper = Point(upper.m_XPos + m_gridUnit, upper.m_YPos + m_gridUnit);
                else
                    upper = Point(upper.m_XPos + m_gridUnit, upper.m_YPos);

                preInterPoint = interPoint;
            }
        }
    } else {
        int x = int((startPoint.m_XPos - m_MinPoint.m_XPos) / m_gridUnit) + 1;
        int y = int((startPoint.m_YPos - m_MinPoint.m_YPos) / m_gridUnit);

        if (fabs((int(startPoint.m_YPos - m_MinPoint.m_YPos) - (startPoint.m_YPos - m_MinPoint.m_YPos))) < ZERO)
            y -= 1;

        Point lower(x * m_gridUnit + m_MinPoint.m_XPos, y * m_gridUnit + m_MinPoint.m_YPos);
        Point preInterPoint = startPoint;

        while (preInterPoint.m_XPos < endPoint.m_XPos) {
            double tempSlope = getSlope(preInterPoint, lower);

            if (slope > tempSlope) {
                double tempY = slope * (lower.m_XPos - preInterPoint.m_XPos) + preInterPoint.m_YPos;
                Point interPoint(lower.m_XPos, tempY);
                double midstX = (preInterPoint.m_XPos + interPoint.m_XPos) / 2;
                double midstY = (preInterPoint.m_YPos + interPoint.m_YPos) / 2;

                int gridX = int((midstX - m_MinPoint.m_XPos) / m_gridUnit) + 1;
                int gridY = int((midstY - m_MinPoint.m_YPos) / m_gridUnit) + 1;

                m_Grids[gridX].m_Grids[gridY].insertSubRoad(gridX, gridY, subRoad);

                double lowerY = slope * (lower.m_XPos - startPoint.m_XPos) + startPoint.m_YPos;

                if (fabs(lowerY - lower.m_YPos) < ZERO)
                    lower = Point(lower.m_XPos + m_gridUnit, lower.m_YPos - m_gridUnit);
                else
                    lower = Point(lower.m_XPos + m_gridUnit, lower.m_YPos);

                preInterPoint = interPoint;

            } else {
                double tempX = (lower.m_YPos - preInterPoint.m_YPos) / slope + preInterPoint.m_XPos;
                Point interPoint(tempX, lower.m_YPos);
                double midstX = (preInterPoint.m_XPos + interPoint.m_XPos) / 2;
                double midstY = (preInterPoint.m_YPos + interPoint.m_YPos) / 2;

                int gridX = int((midstX - m_MinPoint.m_XPos) / m_gridUnit) + 1;
                int gridY = int((midstY - m_MinPoint.m_YPos) / m_gridUnit) + 1;

                m_Grids[gridX].m_Grids[gridY].insertSubRoad(gridX, gridY, subRoad);


                double lowerY = slope * (lower.m_XPos - startPoint.m_XPos) + startPoint.m_YPos;

                if (fabs(lowerY - lower.m_YPos) < ZERO)
                    lower = Point(lower.m_XPos + m_gridUnit, lower.m_YPos - m_gridUnit);
                else
                    lower = Point(lower.m_XPos, lower.m_YPos - m_gridUnit);

                preInterPoint = interPoint;
            }
        }

    }

    /*
        int endX = ( ( endPoint.m_XPos - m_MinPoint.m_XPos ) / m_gridUnit ) + 1;
        int endY = ( ( endPoint.m_YPos - m_MinPoint.m_YPos ) / m_gridUnit ) + 1;
        m_Grids[endX].m_Grids[endY].insertSubRoad(endX,endY, subRoad );

    */

}

/** @brief getSlope

    @todo: document this function
*/
double ComAlgorithm::getSlope(const Point &startPoint, const Point &endPoint) {
    if (fabs(startPoint.m_XPos - endPoint.m_XPos) < ZERO)
        return MAXNUM;

    return (startPoint.m_YPos - endPoint.m_YPos) / (startPoint.m_XPos - endPoint.m_XPos);
}

/** @brief getNearestPoint

    @todo: document this function
*/
void ComAlgorithm::getNearestPoint(const string &strFileName, Roads *const &roads,
                                   Vertexs *const &vertexs, DataPoints *const &dataPoints, INE &ine) {
    if (!strFileName.compare("")) {
        cout << "数据点文件错误！" << endl;
        return;
    }

    //statisticSubRoadofGtids ( roads, vertexs, dataPoints );

    ifstream ifile(strFileName);
    string strBuffer;

    ofstream ofile("queryResult.txt");
    ofstream roadfile("subRoad.txt");

    //int lineNum = -1;

    while (getline(ifile, strBuffer)) {
        if (strBuffer.empty())
            continue;

        vector<string> buffers;
        boost::split(buffers, strBuffer, boost::is_any_of(" "));
        Point point(stod(buffers[1]), stod(buffers[2]));

        //lineNum = stoi ( buffers[0] );

        int gridX = int((point.m_XPos - m_MinPoint.m_XPos) / m_gridUnit) + 1;
        int gridY = int((point.m_YPos - m_MinPoint.m_YPos) / m_gridUnit) + 1;

        vector<SubRoad> subRoads = m_Grids[gridX].m_Grids[gridY].subRoads;

        if (subRoads.size() == 0) {
            ofile << "查询点(" << point.m_XPos << "," << point.m_YPos << ")不在任何道路上" << endl;
            continue;
        }

        Point nearestPoint;
        double vDistance = MAXNUM;
        SubRoad curSubRoad;

        for (SubRoad subRoad : subRoads) {
            roadfile << "road ID:" << subRoad.m_Rid << "\tnode code:" << subRoad.m_Oid << "\t";
            double dx = subRoad.m_End.m_XPos - subRoad.m_Start.m_XPos;
            double dy = subRoad.m_End.m_YPos - subRoad.m_Start.m_YPos;
            double b = dx * subRoad.m_Start.m_YPos - dy * subRoad.m_Start.m_XPos;
            //double vDist = fabs((dy * point.m_XPos - dx * point.m_YPos + dx * subRoad.m_Start.m_YPos - \
                                 dy * subRoad.m_Start.m_XPos) / sqrt(dx * dx + dy * dy));
            double vDist = dy * point.m_XPos - dx * point.m_YPos + b;
            if (fabs(vDistance) > fabs(vDist)) {
                vDistance = vDist;
                curSubRoad = subRoad;

            }
        }

        if (curSubRoad.m_IsVertex) {
            roadfile << "vertex grid coordinate:(" << curSubRoad.m_Oid / MAXLEN << "," << curSubRoad.m_Oid % MAXLEN
                     << ")" << "\n";
            Entity entity = ine.searchNNOfQueryPoint(INEPoint(point.m_XPos, point.m_YPos));
            nearestPoint = Point(entity.m_X, entity.m_Y);
        } else {
            roadfile << "node grid coordinate:(" << curSubRoad.m_Oid / MAXLEN << "," << curSubRoad.m_Oid % MAXLEN
                     << ")" << "\n";
            DataPoint dp = dataPoints[curSubRoad.m_Oid / MAXLEN].m_DataPoints[curSubRoad.m_Oid % MAXLEN];
            nearestPoint = dp.m_Point;

        }

        ofile << "The nearest neighbor of query point " << buffers[0] << " (" << point.m_XPos << "," << point.m_YPos
              << ") is (" << nearestPoint.m_XPos <<
              "," << nearestPoint.m_YPos << ")" << endl;
    }

    ifile.close();
    ofile.close();
}

double ComAlgorithm::getDistanceOfTwoPoint(const Point &startPoint, const Point &endPoint) {
    return sqrt(pow(startPoint.m_XPos - endPoint.m_XPos, 2) + pow(startPoint.m_YPos - endPoint.m_YPos, 2));
}

string &ComAlgorithm::trimString(string &s) {
    if (s.empty())
        return s;

    s.erase(0, s.find_first_not_of("\t"));
    s.erase(s.find_last_not_of("\t") + 1);

    s.erase(0, s.find_first_not_of("\r"));
    s.erase(s.find_last_not_of("\r") + 1);

    s.erase(0, s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);

    return s;
}

void ComAlgorithm::createQueryData(Roads *const &roads) {
    ofstream ofile("TestData/QueryNode.txt", ios::out | ios::trunc);
    int pointNum = 0;
    bool isEnd = false;

    for (size_t row = 0; row < MAXLEN; row++) {
        if (isEnd == true)
            break;

        for (size_t col = 0; col < MAXLEN; col++) {
            Road road;

            if (roads != nullptr && roads[row].m_Roads != nullptr)
                road = roads[row].m_Roads[col];
            else {
                isEnd = true;
                break;
            }
            if (road.m_Rid == -1) {
                isEnd = true;
                break;
            }

            double distance = road.m_Distance;

            Point point = calPointFromDistance(distance * 0.4, road);
            ofile << pointNum << " " << point.m_XPos << " " << point.m_YPos << "\n";
            pointNum++;


        }
    }
    ofile.close();
}

void ComAlgorithm::NaiveNN(const string &strFileName, Roads *const &roads, DataPoints *const &datapoints) {
    if (strFileName.empty()) {
        cout << "查询点文件不存在！" << endl;
        return;
    }

    ifstream iFile(strFileName);
    ofstream oFile("NaiveQueryResult.data", ios::out | ios::trunc);
    string strBuffer;
    int queryNum = 0;

    while (getline(iFile, strBuffer)) {
        strBuffer = trimString(strBuffer);

        if (strBuffer.empty())
            continue;

        vector<string> buffers;
        boost::split(buffers, strBuffer, boost::is_any_of("\t"));
        Point point(stod(buffers[1]), stod(buffers[2]));
        Point nearestPoint;
        int roadID = -1;
        bool isEnd = false;

        for (size_t row = 0; row < MAXLEN; row++) {
            if (isEnd == true)
                break;

            for (size_t col = 0; col < MAXLEN; col++) {
                if (roads[row].m_Roads != nullptr) {
                    Road road = roads[row].m_Roads[col];

                    if (road.m_Rid == -1) {
                        isEnd = true;
                        break;
                    }

                    Point startPoint = road.m_Begin.m_Point;
                    double y = road.m_Slope * (point.m_XPos - startPoint.m_XPos) + startPoint.m_YPos;
                    double diff = -1;
                    y > point.m_YPos ? diff = y - point.m_YPos : diff = point.m_YPos - y;

                    if (fabs(diff) < 1e-4) {
                        roadID = road.m_Rid + 1;
                        double minLen = MAXNUM;

                        if (road.m_Datas.size() > 0) {
                            for (DataPoint dp : road.m_Datas) {
                                double len = getDistanceOfTwoPoint(dp.m_Point, point);

                                if (minLen > len) {
                                    minLen = len;
                                    nearestPoint = dp.m_Point;

                                }
                            }
                        } else {
                            double startlen = getDistanceOfTwoPoint(road.m_Begin.m_Point, point);
                            double endLen = getDistanceOfTwoPoint(road.m_End.m_Point, point);

                            if (startlen > endLen) {
                                minLen = endLen;
                                nearestPoint = road.m_Begin.m_Point;
                            } else {
                                minLen = startlen;
                                nearestPoint = road.m_End.m_Point;
                            }
                        }
                    }
                }
            }
        }

        oFile << "road ID:" << roadID << " 查询点" << queryNum << ":(" << point.m_XPos << "," << point.m_YPos << ")离点("
              << nearestPoint.m_XPos <<
              "," << nearestPoint.m_YPos << ")最近" << endl;
        queryNum++;
    }

}

int ComAlgorithm::timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {

    if (x->tv_sec > y->tv_sec)
        return -1;

    if ((x->tv_sec == y->tv_sec) && (x->tv_usec > y->tv_usec))
        return -1;

    result->tv_sec = (y->tv_sec - x->tv_sec);
    result->tv_usec = (y->tv_usec - x->tv_usec);

    if (result->tv_usec < 0) {
        result->tv_sec--;
        result->tv_usec += 1000000;
    }

    return result->tv_sec + result->tv_usec / 1000;
}
