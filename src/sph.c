#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <float.h>

/***************/
/* Definitions */
/***************/
#define FAST_ACC_COMP 1

#define MAX(a,b) ( (a)>(b)?(a):(b) )
#define MIN(a,b) ( (a)<(b)?(a):(b) )
#define HYPOT2(x, y) sqrt(sqr(x) + sqr(y))
#define sqr(x) ((x)*(x))
#define POW5(x) ((x)*(x)*(x)*(x)*(x))
#define POW4(x) ((x)*(x)*(x)*(x))

#define TRUE 1
#define FALSE 0
#define POWDIM 4
#define DIM 2

#define CUBICSPLINE 0
#define HIGHERORDER 1

#if CUBICSPLINE == 1
static const double SIGMA = 40.0 / (7.0 * M_PI);
#endif

#if HIGHERORDER == 1
static const double SIGMA = 15309 / (478 * M_PI);
#endif

static const double H = 0.2;            // Cut-off radius

static const double GAMMA = 0.1;        
static const double LAMBDA = 7;

static const double RHO0 = 0.001;
static const double ALPHA = 1.0;
static const double BETA = 2.0;

static const double ETA = 0.005;        // Time step parameter
static const double EPS = 0.01 * 0.2;   // Smoothing length

/**************/
/* Structures */
/**************/
typedef enum{particle, pseudoParticle} nodetype;

typedef struct{
    double m;
    double rho;
    double pressure;
    double x[DIM];
    double v[DIM];
    double a[DIM];
    double a_old[DIM];
    int moved;
    int todelete;
	int id;
} Particle;

typedef struct{
    double lower[DIM];
    double upper[DIM];
} Box;

typedef struct TreeNode{
    Particle p;
    Box box;
    nodetype node;
    struct TreeNode *son[POWDIM];
} TreeNode;

/**************/
/* Prototypes */
/**************/
void outputResults(TreeNode *root, FILE *fp);
void freeTree(TreeNode *root);
void moveParticles(TreeNode *root);
void setFlags(TreeNode *t);
void moveLeaf(TreeNode *t, TreeNode *root);
void repairTree(TreeNode *t);
void timeIntegration(double t, double dt, double Tmax, TreeNode *root, Box box);
void initData(TreeNode **root, Box *domain, int N);
void inputParameters(double *dt, double *Tmax, Box *box, int *N);
void insertTree(Particle *p, TreeNode *t);
void genData(Particle *p, int n);
void updateX(Particle *p, double dt);
void compX(TreeNode *root, double dt);
void compV(TreeNode *root, double dt);
void compP(TreeNode *t, TreeNode *root);
void density(TreeNode *tl, TreeNode *t, int level);
void compDensity(Particle *tl, Particle *t);
void acceleration(TreeNode *tl, TreeNode *t, int level);
void compAcceleration(Particle *tl, Particle *t);
void maxVel(TreeNode *t, double *vmax); 

int outsideBox(TreeNode *t);
int sonNumber(Box *box, Box *sonbox, Particle *p);

double adaptiveTimeStep(TreeNode *t);
double randi();
double urand(double, double);

/********/
/* Main */
/********/
int main(int argc, char *argv[]){
    TreeNode *root;
    Box box;
    double dt;
    double Tmax;
    int N;
	inputParameters(&dt, &Tmax, &box, &N);
    initData(&root, &box, N);
    timeIntegration(0, dt, Tmax, root, box);
    freeTree(root);
    return 0;
}

/*************/
/* Functions */
/*************/

void outputResults(TreeNode *t, FILE *fp){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++)
            outputResults(t->son[i], fp);
        // Operation on *t
        if(t->node == particle){
            fprintf(fp, "%f %f %f %f %f %f %d\n",
			t->p.x[0],t->p.x[1],
			t->p.v[0],t->p.v[1], 
			t->p.pressure, t->p.rho, t->p.id
			);
        }
    }
}

void moveParticles(TreeNode *root){
    setFlags(root);
    moveLeaf(root, root);
    repairTree(root);
}

void setFlags(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            setFlags(t->son[i]);
        }
        // Operation on *t
        t->p.moved = FALSE;
        t->p.todelete = FALSE;
    }
}

// Check if particle is still in the box
int outsideBox(TreeNode *t){
    for(int d=0; d<DIM; d++){
        if((t->p.x[d] < t->box.lower[d]) || (t->p.x[d] > t->box.upper[d])){
            return TRUE; // return TRUE if particle is outside of his box
        }
    }
    return FALSE;
}

// Move particels in tree if they are outside of their boxes
void moveLeaf(TreeNode *t, TreeNode *root){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
                moveLeaf(t->son[i], root);
        }
    // Operation on *t
        if((t->node == particle)&&(t->p.moved == FALSE)){
            t->p.moved = TRUE;
            if(outsideBox(t) == TRUE){
                insertTree(&t->p, root);
                t->p.todelete = TRUE;
            }
        }
    }
}

// Repair tree
void repairTree(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            repairTree(t->son[i]);
        }
        // Operation on *t
        if(t->node != particle){
            int numberofsons = 0;
            int d = -1;
            for(int i=0; i<POWDIM; i++) {
                if(t->son[i] != NULL) {
                    if(t->son[i]->p.todelete == TRUE){
                        free(t->son[i]);
                        t->son[i] = NULL;
                    }
                    else {
                        numberofsons++;
                        d=i;
                    }
                }
            }

            if(numberofsons == 0){
                t->p.todelete = TRUE;
            }else if ((numberofsons == 1) && (t->son[d]->node == particle)){
                t->p = t->son[d]->p;
                t->node = t->son[d]->node;
                free(t->son[d]);
                t->son[d] = NULL;
            }

        }
    }
}

void freeTree(TreeNode *t){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            freeTree(t->son[i]);
        }
        // Operation on *t
        free(t);
    }
}

double distanceDomainBoxParticle2D(TreeNode *td, TreeNode *t){
    double x_min = td->box.lower[0];
    double x_max = td->box.upper[0];
    double y_min = td->box.lower[1];
    double y_max = td->box.upper[1];
    double x = t->p.x[0];
    double y = t->p.x[1];

    if (x < x_min) {
        if (y < y_min)
            return HYPOT2(x_min-x, y_min-y);
        else if (y <= y_max)
            return x_min - x;
        else
            return HYPOT2(x_min-x, y_max-y);
    } else if (x <= x_max) {
        if (y < y_min)
            return y_min - y;
        else if (y <= y_max)
            return 0;
        else
            return y - y_max;
    } else {
        if (y < y_min)
            return HYPOT2(x_max-x, y_min-y);
        else if (y <= y_max)
            return x - x_max;
        else
            return HYPOT2(x_max-x, y_max-y);
    }
}


void updateX(Particle *p, double dt){
    for(int d=0; d<DIM; d++){
        p->x[d] += 0.5 * dt * p->v[d];
    }
};

void updateV(Particle *p, double dt){
    for(int d=0; d<DIM; d++){
        p->v[d] += dt * p->a[d];
    }
};

void compX(TreeNode *t, double dt){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compX(t->son[i], dt);
        }
        // Operationen on *t
        if(t->node == particle){
            updateX(&t->p, dt);
        }
    }
}

void compV(TreeNode *t, double dt){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compV(t->son[i], dt);
        }
        // Operationen on *t
        if(t->node == particle){
            updateV(&t->p, dt);
        }
    }
}

void compRho(TreeNode *t, TreeNode *root){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compRho(t->son[i], root);
        }
        // Operation on *t
        if(t->node == particle){
            t->p.rho = 10e-6;
            density(t, root, 0);
        }
    }
}

void density(TreeNode *tl, TreeNode *t, int level){
    if((t != tl)&&(t != NULL)){
        if( (distanceDomainBoxParticle2D(t, tl) < H) || (level == 0)){
            if(t->node == particle){
                compDensity(&tl->p, &t->p);
            }else{
                for(int i=0; i<POWDIM; i++){
                    density(tl, t->son[i], level+1);
                }
            }
        }
    }
}

void compDensity(Particle *i, Particle *j){
    double r=0;
    for(int d=0; d<DIM; d++){
        r += sqr(j->x[d] - i->x[d]);
    }
	r = sqrt(r);
    
    double u = r/H;
    double W = SIGMA/(H*H);

#if CUBICSPLINE == 1
    if((0 <= u) && (u < 0.5)){
        W *= 6 * u * u * u - 6 * u * u + 1;
    }else{
		if((0.5 <= u) && (u <= 1)){
        	W *= 2 * (1 - u) * (1 - u) * (1 - u);
		}else{
			W *= 0;
		}
    }
#endif

#if HIGHERORDER == 1
	if((0 <= u) && (u < 1/3)){
		W *= POW5(1-u) - 6 * POW5(2/3 - u) + 15 * POW5(1/3 - u);
	}else{
		if((1/3 <= u) && (u < 2/3)){
			W *= POW5(1-u) - 6 * POW5(2/3 - u);
		}else{
			if((2/3 <= u) && (u <= 1)){
				W *= POW5(1-u);
			}else{
				W *= 0;
			}
		}
	}
#endif

    i->rho += j->m * W;
}


void compA(TreeNode *t, TreeNode *root){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compA(t->son[i], root);
        }
        // Operation on *t
        if(t->node == particle){
            for(int d=0; d<DIM; d++){
                t->p.a[d] = 0;
            }
            acceleration(t, root, 0);
        }
    }
}

void acceleration(TreeNode *tl, TreeNode *t, int level){
    if((t != tl)&&(t != NULL)){
        if( (distanceDomainBoxParticle2D(t, tl) < H) || (level == 0)){
            if(t->node == particle){
                compAcceleration(&tl->p, &t->p);
            }else{
                for(int i=0; i<POWDIM; i++){
                    acceleration(tl, t->son[i], level+1);
                }
            }
        }
    }
}

void compAcceleration(Particle *i, Particle *j){
    double r=0;
    for(int d=0; d<DIM; d++){
        r += sqr(i->x[d] - j->x[d]);
    }
	double r2 = r;
	r = sqrt(r);
	double u = r/H;    

#if CUBICSPLINE == 1
    double NABLA_W = 6 * SIGMA/(H*H*H);

    if((0 <= u) && (u < 0.5)){
        NABLA_W *= u * (3 * u -2);
    }else{
		if((0.5 <= u) && (u <= 1)){
	        NABLA_W *= - (1 - u) * (1 - u);
		}else{
			NABLA_W *= 0;
		}
    }
#endif

#if HIGHERORDER == 1
    double NABLA_W = -5 * SIGMA/(H*H*H);
	if((0 <= u) && (u < 1/3)){
		NABLA_W *= POW4(1-u) - 6 * POW4(2/3 - u) + 15 * POW4(1/3 - u);
	}else{
		if((1/3 <= u) && (u < 2/3)){
			NABLA_W *= POW4(1-u) - 6 * POW4(2/3 - u);
		}else{
			if((2/3 <= u) && (u <= 1)){
				NABLA_W *= POW4(1-u);
			}else{
				NABLA_W *= 0;
			}
		}
	}
#endif

	double piij = 0;

	double ci = sqrt( GAMMA * RHO0 * pow(i->rho, GAMMA - 1) );
	double cj = sqrt( GAMMA * RHO0 * pow(j->rho, GAMMA - 1) );
	double cij = 0.5 * (ci + cj);

	double vr = 0;
	for(int d=0; d<DIM; d++){
		vr += (i->v[d] - j->v[d]) * (i->x[d] - j->x[d]);
	}
	double muij = (H * vr) / (r2 + EPS * H * H);

	double rhoij = 0.5 * (i->rho + j->rho);

	if(vr < 0){
		piij = (BETA * muij * muij - ALPHA * cij * muij) / rhoij;
	}else{
		piij = 0;
	}

    for(int d=0; d<DIM; d++){
#if FAST_ACC_COMP == 1
        i->a[d] += -(j->m * ( (j->pressure + i->pressure)/(i->rho * j->rho) + piij) * ((i->x[d] - j->x[d])/r) * NABLA_W);
#endif

#if FAST_ACC_COMP == 0
        i->a[d] += -(j->m * ( ((j->pressure / (pow(i->rho,2-LAMBDA) * pow(j->rho,LAMBDA)) ) + 
							   (i->pressure / (pow(i->rho, LAMBDA) * pow(j->rho, 2-LAMBDA)) )) + piij) *
						 	((i->x[d] - j->x[d])/r) * NABLA_W);
#endif

    }
}

void compP(TreeNode *t, TreeNode *root){
    if(t != NULL){
        for(int i=0; i<POWDIM; i++){
            compP(t->son[i], root);
        }

        // Operation on *t
        if(t->node == particle){
            t->p.pressure = RHO0 * pow(t->p.rho, GAMMA);
        }
    }
}

double adaptiveTimeStep(TreeNode *root){
	double vmax = 0;
	maxVel(root, &vmax);
	double dt = ETA * H / vmax;
	if((dt < 10e-5) || (dt > 0.001)){
		dt = 10e-5;
	}
	return dt;
}

void maxVel(TreeNode *t, double *vmax){
	if(t != NULL){
		for(int i=0; i<POWDIM; i++){
			maxVel(t->son[i], vmax);
		}

		// Operation on *t
		if(t->node == particle){
			for(int d=0; d<DIM; d++){
				if((*vmax) < t->p.v[d]){
					(*vmax) = t->p.v[d];
				}
			}
		}
	}
}

void get_file_name(int k, char* buffer, size_t buflen){
    snprintf(buffer, buflen, "logs/data_%d.log", k);
}

void outputResults2File(TreeNode *root, int n){
    const size_t BUFLEN = 50;
    char file_name[BUFLEN];
    get_file_name(n, file_name, BUFLEN);
    FILE *fp;
    fp = fopen(file_name, "w+");
    outputResults(root, fp);
    fclose(fp);
}

void timeIntegration(double t, double dt, double Tmax, TreeNode *root, Box box){
	
    compRho(root, root);
    compP(root, root);
    compA(root, root);
    compV(root, dt/2);
    compX(root, dt/2);
	moveParticles(root);

    int n=0;
    while (t<Tmax){
		dt = adaptiveTimeStep(root);
        t += dt;

        if(n % 50 == 0){
            outputResults2File(root, n);
        }
		
        compRho(root, root);
        compP(root, root);
        compA(root, root);
        compV(root, dt);
        compX(root, dt);

        moveParticles(root);

        n++;
    }
}

void initData(TreeNode **root, Box *domain, int N){
    Particle *p = (Particle*)malloc(N*sizeof(*p));
    // Generate Data
    genData(p, N);
    
    *root = (TreeNode*)calloc(1, sizeof(TreeNode));
    (*root)->p = p[0];
    (*root)->box = *domain;

    for(int i=1; i<N; i++){
        insertTree(&p[i], *root);
    }
    free(p);
}

void inputParameters(double *dt, double *Tmax, Box *box, int *N){
    *dt = 0.00005;
    *Tmax = 5;
    *N = 10100;
    for(int d=0; d<DIM; d++){
        box->lower[d] = -100;
        box->upper[d] = 100;
    }
}

double randi(){
    return (double)rand() / (double)RAND_MAX ; ///(0,1]
}

double urand(double low, double high){
    return low+randi()*(high-low);
}

void genData(Particle *p, int n){
    time_t t;
    srand((unsigned) time(&t));

    double PARTICLE_DISTANCE = 0.1;
    double eta = 1e-4;

    int N_1 = 100;
    int N_2 = 10000;

    double dx = PARTICLE_DISTANCE;
    double dy = PARTICLE_DISTANCE;

    // Object 1
    double px = 0;
    double py = 0;
    double pos_x = -11.*dx;
    double pos_y = 0.5 * dy * 100 - 0.5 * dy * 10 + 0.5 * dy;

    int c = 0;
    for(int i=0; i<N_1; i++){
        c++;
        p[i].x[0] = pos_x + px;
        p[i].x[1] = pos_y + py;
        px += dx;
        if(c%10 == 0){
            px = 0;
            py += dy;
        }
        p[i].v[0] = 100 + eta*urand(-1,1);
        p[i].v[1] = 0 + eta*urand(-1,1);
        p[i].m = 1;
	    p[i].id = 0;
    }

    // Object 2
    px = 0;
    py = 0;
    pos_x = 0;
    pos_y = 0;
    c = 0;
    int k = 0;
    for(int i=N_1; i<N_1+N_2; i++){
        c++;
        p[i].x[0] = pos_x + px;
        p[i].x[1] = pos_y + py;
        py += dy;
        if(c%100 == 0){
            k++;
            px += dx;
            py = 0;
            if(k%10 == 0){
                px += 5*dx;
            }
        }
        p[i].v[0] = eta*urand(-1,1); 
        p[i].v[1] = eta*urand(-1,1); 
        p[i].m = 1;
	    p[i].id = 1;
    }
}

int sonNumber(Box *box, Box *sonbox, Particle *p){
    int b=0;
    for(int d=DIM-1; d>=0; d--){
        if(p->x[d] < .5*(box->upper[d] + box->lower[d])){
            b = 2*b;
            sonbox->lower[d] = box->lower[d];
            sonbox->upper[d] = .5*(box->upper[d] + box->lower[d]);
        }else{
            b = 2*b+1;
            sonbox->lower[d] = .5*(box->upper[d] + box->lower[d]);
            sonbox->upper[d] = box->upper[d];
        }
    }
    return b;
}

void insertTree(Particle *p, TreeNode *t){
    Box sonbox;
    int b=sonNumber(&t->box, &sonbox, p);

    if(t->son[b] == NULL){
        if(t->node == particle){
            Particle p2 = t->p;
            t->son[b] = (TreeNode*)calloc(1,sizeof(TreeNode));
            t->son[b]->box = sonbox;
            t->son[b]->p = *p;

            t->son[b]->node = particle;
            t->node = pseudoParticle;

            insertTree(&p2, t);
        }else{
            t->son[b] = (TreeNode*)calloc(1,sizeof(TreeNode));
            t->son[b]->box = sonbox;
            t->son[b]->p = *p;

            t->son[b]->node = particle;
        }
    }else{
        insertTree(p, t->son[b]);
    }
}
