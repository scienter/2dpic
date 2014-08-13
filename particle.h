

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _ptclList  {
    float x; 
    float oldX;   
    float y; 
    float oldY;   
    float p1;    //momentum  
    float p2;
    float p3;
    float E1;    
    float B1;    
    float Pr;    
    float Pl;    
    float Sr;    
    float Sl;    
    float index; 
    struct _ptclList *next;
} ptclList;

