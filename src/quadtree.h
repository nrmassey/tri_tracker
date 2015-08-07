/******************************************************************************
** Program : quadtree.h
** Author  : Neil Massey
** Date    : 07/07/09
** Purpose : class to contain a quadtree to represent the hierarchical 
**           triangular mesh
******************************************************************************/

#ifndef QT_H
#define QT_H

#include <stddef.h>
#include <assert.h>
#include <iostream>
#include <algorithm>

template <class T>
class qt_node
{
    friend std::ostream& operator<<(std::ostream& os, qt_node<T>* node)
    {
        os << node->get_data();
        for (int i=0; i<4; i++)
        {
            if (node->get_child(i) != NULL)
                os << node->get_child(i);
        }
        return os;
    }

    public:
    
        /*********************************************************************/ 
            
        qt_node(void)
        {
            for (int i=0; i<4; i++)
                child[i] = NULL;
            parent = NULL;
            nc = 0;
        }

        /*********************************************************************/     
        
        qt_node(qt_node* pparent)
        {
            for (int i=0; i<4; i++)
                child[i] = NULL;
            parent = pparent;
            nc = 0;
        }

        /*********************************************************************/ 
            
        ~qt_node(void)
        {
            for (int i=0; i<4; i++)
                delete child[i];
        }

        /*********************************************************************/
        
        qt_node* add_child(void)
        {
            assert(nc < 4);
            child[nc] = new qt_node(this);
            nc++;
            return child[nc-1];
        }

        /*********************************************************************/
        
        qt_node* add_child(T child_data)
        {
            assert(nc < 4);     // can't add any more children after array is full
            child[nc] = new qt_node(this);
            child[nc]->data = child_data;
            nc++;
            return child[nc-1];
        }

        /*********************************************************************/
        
        qt_node<T>* get_children(void)
        {
            return child;
        }

        /*********************************************************************/
        
        qt_node<T>* get_child(int child_number)
        {
            if (child_number < 0 || child_number > 3)
                return NULL;
            else
                return child[child_number];
        }

        /*********************************************************************/
        
        qt_node<T>* get_parent(void)
        {
            return parent;
        }

        /*********************************************************************/
        
        T* get_data(void)
        {
            return &data;
        }
        
        /*********************************************************************/
    
        bool is_leaf(void)
        {
            bool is_leaf = true;
            for (int i=0; i<4; i++)
                is_leaf = is_leaf && (child[i] == NULL);
            return is_leaf;
        }

        /*********************************************************************/
        
        int get_level(void)
        {
            // get the depth the node occurs at
            int l = 0;
            qt_node<T>* cn = parent;
            while (cn != NULL)
            {
                l++;
                cn = cn->parent;
            }
            return l;
        }
        
        /*********************************************************************/
        
        int get_max_level(qt_node<T>* node)
        {
            // get the maximum number of levels under this node
            if (node == NULL)
                return 0;
            else
            {
                std::vector<int> mD(4);
                for (int i=0; i<4; i++)
                    mD[i] = get_max_level(node->get_child(i));
                std::sort(mD.begin(), mD.end());
                return mD[3]+1;
            }           
        }
        
    /*************************************************************************/

    private:
        T data;
        qt_node<T>* parent;
        qt_node<T>* child[4];
        int nc;                 // current number of children
};

/*****************************************************************************/

template <class T>
class quadtree
{
    friend std::ostream& operator<<(std::ostream& os, quadtree<T>& QT)
    {
        os << QT.get_root();
        return os;
    }

    public:

        /*********************************************************************/ 

        quadtree(void) : root(NULL)
        {
            root = new qt_node<T>;
        }

        /*********************************************************************/

        quadtree(T root_data) : root(NULL)
        {
            root = new qt_node<T>;
            *(root->get_data()) = root_data;
        }

        /*********************************************************************/

        ~quadtree(void)
        {
            delete root;
        }       

        /*********************************************************************/

        qt_node<T>* get_root(void)
        {
            return root;
        }
        
        /*********************************************************************/

        int get_max_level(void)
        {
            return root->get_max_level(root);
        }
        
        /*********************************************************************/

        std::list<qt_node<T>* > get_all_nodes_at_level(int level)
        {
            std::list<qt_node<T>* > node_list;
            get_nodes_at_level(root, level, &node_list);
            return node_list;
        }

        /*********************************************************************/
        
    private:
    
        /*********************************************************************/
                
        void get_nodes_at_level(qt_node<T>* current_node, int level,
                                std::list<qt_node<T>* >* node_list)
        {
            // check whether the current node equals the level
            if (current_node->get_level() == level)
                node_list->push_back(current_node);

            // recursively get the child nodes
            for (int i=0; i<4; i++)
                if (current_node->get_child(i) != NULL)
                    get_nodes_at_level(current_node->get_child(i), level, node_list);
        }
        
        /*********************************************************************/
    
        qt_node<T>* root;
};

#endif
