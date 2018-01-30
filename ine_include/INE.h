#ifndef INE_H
#define INE_H

#include "INEClassType.h"
#include "CodeBeauty.h"

class INE
{
  public:
    INE();
    virtual ~INE();

  public:
    void readVertexFile(string strFileName);
    void readRoadAndEntity(string strFileName);
    Entity searchNNOfQueryPoint(const INEPoint& queryPoint);


  protected:
    Entity getEntity(const INERoad& road,const double distance);
    INERoad searchRoadCoverQueryPoint(const INEPoint& queryPoint);
    double calDistanceOf2Point(const INEPoint& point1,const INEPoint& point2);

  private:
    map<int,INEVertex> m_Vertexs;
    map<int,Entity> m_Entitys;
    map<int,INERoad> m_Roads;

    //VV m_VVs[MAX_NUM];

    map<int,map<int,set<Data>>> m_VVs;

};

#endif // INE_H
