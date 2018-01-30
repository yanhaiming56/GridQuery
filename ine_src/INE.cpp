#include "../ine_include/INE.h"
#include "../ine_include/INEClassType.h"


INE::INE() {
    //ctor
}

INE::~INE() {
    //dtor
}

void INE::readVertexFile(string strFileName) {
    ifstream ifile(strFileName);
    string strLine;
    while (getline(ifile, strLine)) {
        vector<string> strBuffers;
        boost::split(strBuffers, strLine, boost::is_any_of(" "));
        INEVertex v(stoi(strBuffers[0]), stod(strBuffers[1]), stod(strBuffers[2]));
        m_Vertexs[stoi(strBuffers[0])] = v;
        m_VVs[stoi(strBuffers[0])] = map<int, set<Data>>();
    }

    ifile.close();
}

void INE::readRoadAndEntity(string strFileName) {
    ifstream ifile(strFileName);
    string strLine;
    int64_t road_code = 0;
    int64_t entity_code = 0;
    bool isRoad = true;
    INERoad *curRoad = nullptr;

    while (getline(ifile, strLine)) {
        vector<string> strBuffers;
        boost::split(strBuffers, strLine, boost::is_any_of(" "));

        if (isRoad == true) {
            m_Roads[road_code] = INERoad(road_code, stoi(strBuffers[0]), stoi(strBuffers[1]), stod(strBuffers[2]));
            curRoad = &m_Roads[road_code];
            m_Vertexs[curRoad->m_StartId].m_Roads.insert(curRoad->m_Rid);
            m_Vertexs[curRoad->m_EndId].m_Roads.insert(curRoad->m_Rid);

            map<int, set<Data>>::iterator iter_datas;
            map<int, set<Data>> map_data = m_VVs[curRoad->m_StartId];
            iter_datas = map_data.find(curRoad->m_EndId);
            if (iter_datas == map_data.end())
                m_VVs[curRoad->m_StartId][curRoad->m_EndId] = set<Data>();

            map_data = m_VVs[curRoad->m_EndId];
            iter_datas = map_data.find(curRoad->m_StartId);
            if (iter_datas == map_data.end())
                m_VVs[curRoad->m_EndId][curRoad->m_StartId] = set<Data>();

            if (stoi(strBuffers[3]) != 0)
                isRoad = false;

            road_code++;
        } else {
            for (size_t i = 1; i < strBuffers.size(); i += 2) {
                Entity entity = getEntity(*curRoad, stod(strBuffers[i]));
                m_Entitys[entity_code] = getEntity(*curRoad, stod(strBuffers[i]));
                m_Entitys[entity_code].m_Oid = entity_code;
                entity = m_Entitys[entity_code];
                curRoad->m_Entitys.insert(entity_code);
                m_VVs[curRoad->m_StartId][curRoad->m_EndId].insert(
                        Data(entity_code, stod(strBuffers[i])));
                m_VVs[curRoad->m_EndId][curRoad->m_StartId].insert(
                        Data(entity_code, curRoad->m_Distance - stod(strBuffers[i])));
                entity_code++;
            }
            isRoad = true;
        }
    }
    ifile.close();
}

Entity INE::getEntity(const INERoad &road, const double distance) {
    Entity entity;
    INEVertex start_v = m_Vertexs[road.m_StartId];
    INEVertex end_v = m_Vertexs[road.m_EndId];
    double dx = end_v.m_X - start_v.m_X;
    double dy = end_v.m_Y - start_v.m_Y;

    //垂直线
    if (fabs(dx) < ZERO) {
        end_v.m_Y > start_v.m_Y ? entity.m_Y = start_v.m_Y + distance : entity.m_Y = start_v.m_Y - distance;
        entity.m_X = start_v.m_X;
        return entity;
    }

    //水平线
    if (fabs(dy) < ZERO) {
        end_v.m_X > start_v.m_X ? entity.m_X = start_v.m_X + distance : entity.m_X = start_v.m_X - distance;
        entity.m_Y = start_v.m_Y;
        return entity;
    }

    double b = dx * start_v.m_Y - dy * start_v.m_X;
    double x = start_v.m_X - fabs(dx * distance) / sqrt(dx * dx + dy * dy);

    if (fabs(x - start_v.m_X) + fabs(x - end_v.m_X) > fabs(dx))
        x = start_v.m_X + fabs(dx * distance) / sqrt(dx * dx + dy * dy);

    double y = (dy * x + b) / dx;
    entity.m_X = x;
    entity.m_Y = y;
    return entity;
}

INERoad INE::searchRoadCoverQueryPoint(const INEPoint &queryPoint){
    for (map<int, INERoad>::iterator iter = m_Roads.begin(); iter != m_Roads.end(); iter++) {
        INERoad road = iter->second;
        INEVertex start_v = m_Vertexs[road.m_StartId];
        INEVertex end_v = m_Vertexs[road.m_EndId];
        double dx = end_v.m_X - start_v.m_X;
        double dy = end_v.m_Y - start_v.m_Y;
        double b = dx * start_v.m_Y - dy * start_v.m_X;
        double result = dy * queryPoint.m_X - dx * queryPoint.m_Y + b;

        if (fabs(result) < ZERO)
            return road;
    }
}

Entity INE::searchNNOfQueryPoint(const INEPoint &queryPoint) {
    priority_queue<VisitedNode> pq;
    map<int, bool> visitedRoad;

    INERoad cover_road = searchRoadCoverQueryPoint(queryPoint); //get the road coverying the query point
    set<Data> datas = m_VVs[cover_road.m_StartId][cover_road.m_EndId];
    //the road of coverying the query point contains entity
    if (datas.size() > 0) {
        double minDist = MAX_VAL;
        int betterID = -1;
        for (set<Data>::iterator iter = datas.begin(); iter != datas.end(); iter++) {
            double curDist = calDistanceOf2Point((INEPoint) m_Entitys[iter->m_Oid], queryPoint);
            if (minDist > curDist) {
                minDist = curDist;
                betterID = iter->m_Oid;
            }

        }
        return m_Entitys[betterID];

    }

    visitedRoad[cover_road.m_Rid] = true;
    INEVertex start_v = m_Vertexs[cover_road.m_StartId];
    INEVertex end_v = m_Vertexs[cover_road.m_EndId];
    pq.push(VisitedNode(start_v.m_Vid, calDistanceOf2Point((INEPoint) start_v, queryPoint)));
    pq.push(VisitedNode(end_v.m_Vid, calDistanceOf2Point((INEPoint) end_v, queryPoint)));


    set<VisitedEntity> ves; //record the entity which has contained by visited Road
    double entityDist = MAX_VAL; //record the maximum distance about the visited entity from query point

    VisitedNode firstNode = pq.top();
    pq.pop();
    while (firstNode.m_Distance < entityDist) {
        map<int, set<Data>> vds = m_VVs[firstNode.m_Vid]; //get end vertex
        for (map<int, set<Data>>::iterator iter = vds.begin(); iter != vds.end(); iter++) {
            start_v = m_Vertexs[firstNode.m_Vid];
            end_v = m_Vertexs[iter->first];
            int visitedRID = -1;
            bool isEnd = false;
            for (int s_rid : start_v.m_Roads) {
                if (isEnd == true)
                    break;
                for (int e_rid:end_v.m_Roads)
                    if (s_rid == e_rid) {
                        visitedRID = s_rid;
                        isEnd = true;

                    }
            }

            if (visitedRoad.find(visitedRID) != visitedRoad.end())
                continue;
            visitedRoad[visitedRID] = true; //flag all visited road from fist-node to all adjcent node

            set<Data> dts = iter->second; //get the entity about min distance from start vertex
            if (dts.size() > 0) {
                set<Data>::iterator iter = dts.begin(); //only get the nearest entity from vertex
                ves.insert(VisitedEntity(iter->m_Oid, firstNode.m_Distance + iter->m_Distance));
                entityDist = ves.rbegin()->m_Distance;
            }

            INERoad curRoad = m_Roads[visitedRID];

            VisitedNode tempVN(iter->first, firstNode.m_Distance + curRoad.m_Distance);
            pq.push(tempVN);
        }
        firstNode = pq.top();
        pq.pop();
    }
    return m_Entitys[ves.begin()->m_Oid];


}

double INE::calDistanceOf2Point(const INEPoint &point1, const INEPoint &point2) {
    double dx = point1.m_X - point2.m_X;
    double dy = point1.m_Y - point2.m_Y;
    return sqrt(dx * dx + dy * dy);
}



