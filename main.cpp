#include <iostream>

#include "ComAlgorithm.h"
#include "ine_include/INE.h"

using namespace std;
//0.01 TestData/CalNode.txt TestData/CalRoadData.txt

//int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

int main(int argv, char **argc) {
    INE ine;
    timeval startTime, endTime, diffTime;
    ComAlgorithm CA;
    string strFileName;
    ofstream ofile("executTime.data");

    double iGridUnit = 0;
    cout << "Pleas input the unit length of grid!" << endl;
    cin >> iGridUnit;
    CA.setGridUnit(iGridUnit);

    //获取顶点数据
    strFileName.clear();
    cout << "Please input the vertex file" << endl;
    strFileName = "TestData/CalNode.txt";
    Vertexs *vertexs = nullptr;
    gettimeofday(&startTime, 0);
    CA.getVertexs(strFileName, vertexs);
    ine.readVertexFile(strFileName);
    gettimeofday(&endTime, 0);
    CA.timeval_subtract(&diffTime, &startTime, &endTime);
    ofile << "构建顶点索引耗时：" << diffTime.tv_sec * 1000 + diffTime.tv_usec / 1000 << "ms" << endl;

    //获取道路和道路数据点信息
    strFileName.clear();
    cout << "Please input road and data relation file." << endl;
    strFileName = "TestData/CalRoadData.txt";
    Roads *roads = nullptr;
    DataPoints *dataPoints = nullptr;
    gettimeofday(&startTime, 0);
    CA.getRoadOfPoints(strFileName, dataPoints, roads, vertexs);
    ine.readRoadAndEntity(strFileName);
    gettimeofday(&endTime, 0);
    CA.timeval_subtract(&diffTime, &startTime, &endTime);
    ofile << "构建道路及数据点索引耗时：" << diffTime.tv_sec * 1000 + diffTime.tv_usec / 1000 << "ms" << endl;

    //CA.createQueryData(roads);

    CA.statisticSubRoadofGtids(roads, vertexs, dataPoints);

    cout << "<-------Grid查询开始----------->" << endl;
    strFileName.clear();
    //strFileName = "TestData/QueryNode1.data";
    strFileName = "TestData/QueryNode.txt";
    gettimeofday(&startTime, 0);
    CA.getNearestPoint(strFileName, roads, vertexs, dataPoints,ine);
    gettimeofday(&endTime, 0);
    CA.timeval_subtract(&diffTime, &startTime, &endTime);
    ofile << "计算Grid查询点最近邻耗时：" << diffTime.tv_sec * 1000 + diffTime.tv_usec / 1000 << "ms" << endl;
    cout << "计算Grid查询点最近邻耗时：" << diffTime.tv_sec * 1000 + diffTime.tv_usec / 1000 << "ms" << endl;
    cout << "<-------Grid查询结束----------->" << endl;


    if (vertexs != nullptr) {
        for (size_t i = 0; i < MAXLEN; i++)
            delete[] vertexs[i].m_Vertexs;
        delete[] vertexs;
    }
    if (dataPoints != nullptr) {
        for (size_t i = 0; i < MAXLEN; i++)
            delete[] dataPoints[i].m_DataPoints;
        delete[] dataPoints;
    }
    if (roads != nullptr) {
        for (size_t i = 0; i < MAXLEN; i++)
            delete[] roads[i].m_Roads;
        delete[] roads;
    }
    ofile.close();

    return 0;
}


/*int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {

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
}*/
