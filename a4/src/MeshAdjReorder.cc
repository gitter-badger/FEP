#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "MeshAdjReorder.h"

bool hasNode(apf::Mesh* m, apf::MeshEntity* e)
{
    return m->getShape()->countNodesOn(m->getType(e)) > 0;
}

struct sortVertsFunctor {
    apf::Mesh* m;
    bool operator()(apf::MeshEntity *lhs, apf::MeshEntity *rhs) {
    apf::Adjacent adj;
    this->m->getAdjacent(lhs, 1, adj);
    int lhs_count = adj.getSize();
    this->m->getAdjacent(rhs, 1, adj);
    int rhs_count = adj.getSize();
    return (lhs_count < rhs_count);
    }
};

void adjReorder(apf::Mesh* m, 
            apf::FieldShape* shape, 
            int components,
            apf::Numbering* & nodeNums,
            apf::Numbering* & faceNums) 
{
    //std::cout << "inside" << std::endl;
    //std::cout << nodeNums << std::endl;
    //std::cout << faceNums << std::endl;
    //std::cout << m << std::endl;
    //std::cout << shape << std::endl;
    //std::cout << "inside" << std::endl;

    apf::MeshEntity* e;
    apf::MeshIterator* it;
    int node_label = 0;
    int dim = m->getDimension();

    for(int ii = 0; ii <= dim; ++ii){
        if(shape->hasNodesIn(ii)){
            it = m->begin(ii);
            while((e = m->iterate(it))) {
                std::cout << e << std::endl;
                if(m->isOwned(e)) {
                    node_label += shape->countNodesOn(m->getType(e)) * components;
                }
            }
        }
    }
    //std::cout << node_label << std::endl;
    //std::cout << "Mesh: " << m << std::endl;
    //std::cout << "Shape: " << shape << std::endl;

  // deduct any fixed entities from the locally-owned total
    
    unsigned int min_order = 400;
    apf::MeshEntity* start_vert = NULL;
    it = m->begin(0);
    while((e = m->iterate(it))){
        std::cout << "inside" << std::endl;
        apf::ModelEntity* model_ent = m->toModel(e);
        int dim = m->getModelType(model_ent);
        std::cout << dim << std::endl;
        if(dim != 0) {
            std::cout << "continue" << std::endl;
            continue;
        }
        apf::Adjacent adj;
        m->getAdjacent(e, 1, adj);
        //this keeps the first lowest order is sees
        std::cout << adj.getSize() << std::endl;
        if(adj.getSize() < min_order) {
            min_order =  adj.getSize();
            start_vert = e; 
            std::cout << "start" << start_vert << std::endl;
        }
    }
    //begin the node relabeling by creating an empty ordering
    //std::cout << "start_vert: " << start_vert << std::endl;
    //std::cout << "Mesh: " << m << std::endl;
    //std::cout << "Shape: " << shape << std::endl;
    //return;
    std::queue<apf::MeshEntity*> q;
    q.push(start_vert);
    std::set<apf::MeshEntity*> in_q;
    in_q.insert(start_vert);
    std::set<apf::MeshEntity*>::iterator set_it;
    std::vector<apf::MeshEntity*> ent_list;

    signed int face_label = m->count(2);
    //we stop at zero because we labeled every node
    //everything left in the queue will already have been labeled
    while( q.size() > 0 && node_label != 0) {
        //std::cout << "Mesh: " << m << std::endl;
        //std::cout << "Shape: " << shape << std::endl;
        apf::MeshEntity* current_entity;
        current_entity = q.front();
        q.pop();
        set_it = in_q.find(current_entity);
        if(set_it != in_q.end()) {
            in_q.erase(set_it);
            //std::cout << "erased element from set" << std::endl;
        }
        //std::cout << "inside" << std::endl;

        //assuming that there is only one dof per entity
        if(!apf::isNumbered(nodeNums, current_entity,0,0)) {
            //std::cout << "nodeNUms " << nodeNums << std::endl;
            std::cout << "entity " << current_entity << std::endl;
            //std::cout << node_label << std::endl;
            apf::number(nodeNums, current_entity, 0, 0, node_label);
            //std::cout << "labeled entity as: " << node_label << std::endl;
            --node_label;
        } else {
            continue;
        }

        //now we want to add any unlabled adjacent mesh entities to the queue
        //additions to queue will only be handled  by looking at adjacencies from
        //vertices
        if(apf::getDimension(m, current_entity) == 0) {
            apf::MeshEntity* vert = current_entity;
            apf::Adjacent edge_adj;
            m->getAdjacent(vert, 1, edge_adj);
            int edge_adj_size = edge_adj.getSize();
            for(int ii = 0; ii < edge_adj_size; ++ii) {
        apf::MeshEntity* curr_edge = edge_adj[ii];
        //get the face adjacency for this edge and see if it is labeled
        apf::Adjacent face_adj;
        m->getAdjacent(edge_adj[ii], 2, face_adj);
        int face_adj_size = face_adj.getSize();
        for(int jj = 0; jj < face_adj_size; ++jj) {
            if(!apf::isNumbered(faceNums, face_adj[jj], 0, 0)) {
            apf::number(faceNums, face_adj[jj], 0 ,0, face_label);
            --face_label;
            }
            //check if the face has a node to label, then queue it if exists
            if(hasNode(m, face_adj[jj])) { // we will never have a dof on a face
            //queue the face node
            //std::cout << "trying to queue a face node" << std::endl;
            }
        }
        //get the other end of the edge
        apf::MeshEntity* other_vert = apf::getEdgeVertOppositeVert(m, curr_edge, vert);
        if(hasNode(m, curr_edge)) {
            //check that the other_vert is labeled or in queue and edge node not labeled
            set_it = in_q.find(other_vert);
            bool is_in_q = (set_it != in_q.end());
            if( (apf::isNumbered(nodeNums, other_vert, 0, 0) || is_in_q ) &&
                !apf::isNumbered(nodeNums, curr_edge, 0, 0) ) {
            apf::number(nodeNums, curr_edge, 0 , 0, node_label);
            --node_label;
            } else {
            //check that this edge is not labeled or in queue
            set_it = in_q.find(curr_edge);
            if(!apf::isNumbered(nodeNums, curr_edge, 0, 0) && set_it == in_q.end()) {
                q.push(curr_edge);
                in_q.insert(curr_edge);
            }
            //add the other vertex to the list
            ent_list.push_back(other_vert);
            }
        } else {
            if(!apf::isNumbered(nodeNums, other_vert, 0, 0)) {
            ent_list.push_back(other_vert);
            }
        }
            }// finished loop over all edges coming into vertex
            //add all of the elements in the list to the queue
            if(ent_list.size() > 0) {
        //use the functor to create a closure over the mesh for ranking the vector
        struct sortVertsFunctor functor;
        functor.m = m;
        std::sort(ent_list.begin(), ent_list.end(), functor);
        //add each element of the vector to the queue
        //print out the ent vector and the queue
        for(std::vector<apf::MeshEntity*>::iterator it = ent_list.begin();
            it != ent_list.end(); ++it) {
            //check that it is not already in queue
            set_it = in_q.find(*it);
            if(set_it == in_q.end()) {
            q.push(*it);
            in_q.insert(*it);
            } else {
            ////std::cout << "already in queue" << std::endl;
            }
        }
        ent_list.clear();
            }
        } else {
            ////std::cout << "not a vertex" << std::endl;
        } // end of vertex checking
        if(q.size() > 500 ) {
            //std::cout << "bad termination" << std::endl;
            break;
        }
        } //end of while loop

        //make sure that we did not leave an unlabeled node in the queue
        while(!q.empty()) {
        assert(apf::isNumbered(nodeNums, q.front(), 0, 0) == 1);
        q.pop();
    }
}

