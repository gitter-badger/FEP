#include <queue>
#include <stdio>
#include <stdlib>

class MeshWrapper {
  MeshWrapper(apf::MeshEntity* mE, int val) 
    {this->me = mE; this->dof = val;}
  ~MeshWrapper();
  operator< (const &lhs, const &rhs) {return (lhs.dof < rhs.dof);}

  apf::MeshEntity *me;
  int dof;
};
do
{
  /* code */
} while (/* condition */);

/*Pseudo Code*/
Renumber {
  /*there is only one number per node*/
  apf::Numbering* nodeNums = apf::createNumbering(mesh, "nodeNums", /
    m->getShape(), 1);
  int label = n+1;
  MeshEntity* first_node_holder = getStart();
  std::queue<MeshEntity*> q;
  q.push(first_node_holder);

  while(!q.empty()){
    label = label - 1;
    MeshEntity* current_entity = q.pop();
    /*We use the pop until we have an unlabled node version
    because implementing a set of labeld nodes is equivalent to
    recreating the array of mesh pointers*/
    /*(-1) currently means unlabled*/
    while(apf::getNumber(nodeNums, current_entity, 0, 0) != -1){
      current_entity = q.pop();
    }
    /*check if the current entity is */
    if(apf::FieldShape::hasNodesIn(1)) {
      if(mesh->getType(current_entity) == apf::Mesh::EDGE) {
        /*get the other edge that is unlabeled, because one will always
        be labeled in this scheme*/
        apf::Downward down;
        int num_down_verts = mesh->getDownward(current_entity, 0, down);
        /*add all unlabeled node to the queue, this should only be one*/
        for(int i = 0; i < num_down_verts; ++i) {
          if(apf::isNumbered(nodeNums, down[i], 0, 0) ) {
            /*pass*/
          } else {
            q.push_back(down[i]);

            std::cout << "added unlabeled node from edge: " << /
            apf::getNumber(nodeNums, current_entity, 0, 0) /
            << std::endl;
          }
        }

      } else if {

      } else {
        std::cout << "something is wrong. we should have identified /
        the node by now" << std::endl;
      }

    } else {
      /*this means we jump over the edges to vertices on the other sides*/
      apf::Adjacent adj;
      /*get the edge around*/
      mesh->getAdjacent(current_entity, 1, adj);
      apf::MeshEntity* adj_edge;
      std::std::vector<MeshWrapper> wrapped;
      for(int i = 0; i < adj.getSize(); ++i) {
        apf::MeshEntity* other_vert = /
          apf::getEdgeVertOppositeVert(mesh, adj_edge, current_entity);
        apf::Adjacent other_adj;
        mesh->getAdjacent(other_vert)
        wrapped.push_back(MeshWrapper(other_vert, ))

      }
      /*sort the vector*/
    }

  }
}