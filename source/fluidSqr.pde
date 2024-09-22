final int N = 200;
final int iter = 16;

int coords(int x, int y){
  x = constrain(x, 0, N - 1);
  y = constrain(y, 0, N-1);
  return x + y*N;
}

class FluidSqr {
  int size;
  float dt;
  float diff;
  float visc;
  
  float[] s;
  float[] density;
  
  float[] Vx;
  float[] Vy;
  
  float[] Vx0;
  float[] Vy0;
  
  FluidSqr (float dt , float diffusion, float viscosity){
    this.size = N;
    this.dt = dt;
    this.diff = diffusion;
    this.visc = viscosity;
    this.s = new float[N*N];
    this.density = new float[N*N];
    this.Vx = new float[N*N];
    this.Vy = new float[N*N];
    this.Vx0 = new float[N*N];
    this.Vy0 = new float[N*N];
  }
  
  void renderD() {
    //colorMode(HSB, 255);
    for (int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        float x = i * SCALE;
        float y = j * SCALE;
        float d = this.density[coords(i,j)];
        //fill((frameCount%255)*d,d,0);
        fill(d*30);
        //fill((d + 50) % 255,200,d);
        //stroke(1);
        noStroke();
        square(x, y, SCALE);
      }
    }
  }
  
  
  void step(){
    float visc     = this.visc;
    float diff     = this.diff;
    float dt       = this.dt;
    float[] Vx      = this.Vx;
    float[] Vy      = this.Vy;
    float[] Vx0     = this.Vx0;
    float[] Vy0     = this.Vy0;
    float[] s       = this.s;
    float[] density = this.density;
    
    
    
    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);
    
    project(Vx0, Vy0, Vx, Vy);
    
    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);
    
    project(Vx, Vy, Vx0, Vy0);
    
    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
}
  
  
  
  
  
  void addDensity(int x, int y, float amount){
    int c = coords(x,y);
    this.density[c] += amount;
  }
  
  void addVelocity(int x, int y, float amountX, float amountY){
    int c = coords(x,y);
    this.Vx[c] += amountX;
    this.Vy[c] += amountY;
  }
  
  void diffuse (int b, float[] x, float[] x0, float diff, float dt){
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 4 * a);
}
  void lin_solve(int b, float[] x, float[] x0, float a, float c){
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
      for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            x[coords(i, j)] =
              (x0[coords(i, j)]
              + a*(x[coords(i+1, j)]
              +x[coords(i-1, j)]
              +x[coords(i  , j+1)]
              +x[coords(i  , j-1)]                    
              )) * cRecip;
                }
            }
        set_bnd(b, x);
    }
  }
  
  void set_bnd(int b, float[] x){
      //switched order of b == 2 and b == 1
        for(int i = 1; i < N - 1; i++) {
            x[coords(i, 0 )] = b == 2 ? -x[coords(i, 1)] : x[coords(i, 1)];
            x[coords(i, N-1)] = b == 2 ? -x[coords(i, N-2)] : x[coords(i, N-2)];
        }
        for(int j = 1; j < N - 1; j++) {
            x[coords(0  , j)] = b == 1 ? -x[coords(1  , j)] : x[coords(1  , j)];
            x[coords(N-1, j)] = b == 1 ? -x[coords(N-2, j)] : x[coords(N-2, j)];
        }
    
    //was 0.33f instead of 0.5
    x[coords(0, 0)] = 0.5f * (x[coords(1, 0)]
                      + x[coords(0, 1)]);
                                  
    x[coords(0, N-1)] = 0.5f * (x[coords(1, N-1)]
                      + x[coords(0, N-2)]);

    x[coords(N-1, 0)] = 0.5f * (x[coords(N-2, 0)]
                      + x[coords(N-1, 1)]);
                                  
    x[coords(N-1, N-1)] = 0.5f * (x[coords(N-2, N-1)]
                      + x[coords(N-1, N-2)]);

  }
  
  void project(float[] velocX, float[] velocY, float[] p, float[] div) {
    for (int j = 1; j < N - 1; j++) {
      for (int i = 1; i < N - 1; i++) {
        div[coords(i, j)] = -0.5f*(
            velocX[coords(i+1, j)]
            -velocX[coords(i-1, j)]
            +velocY[coords(i , j+1)]
            -velocY[coords(i  , j-1)])/N;
            p[coords(i, j)] = 0;
            }
        }
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6); // was 6, i believe its the amount of dimentions
    
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[coords(i, j)] -= 0.5f * (  p[coords(i+1, j)]
                                                -p[coords(i-1, j)]) * N;
                velocY[coords(i, j)] -= 0.5f * (  p[coords(i, j+1)]
                                                -p[coords(i, j-1)]) * N;
            }
        }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
  }
  void advect(int b, float[] d, float[] d0,  float[] velocX, float[] velocY, float dt){
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[coords(i, j)];
                tmp2 = dty * velocY[coords(i, j)];

                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;

                //was 0.5f instead of 1
                if(x < 1f) x = 1f; 
                if(x > Nfloat + 1f) x = Nfloat + 1f; 
                i0 = x; 
                i1 = i0 + 1.0f;
                if(y < 1f) y = 1f; 
                if(y > Nfloat + 1f) y = Nfloat + 1f; 
                j0 = y;
                j1 = j0 + 1.0f; 

                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;

                
              //i0,i1,j0,j1 were turned to int for the next part
              int i0i = (int) i0;
              int i1i = (int) i1;
              int j0i = (int) j0;
              int j1i = (int) j1;
                
                d[coords(i, j)] = 
                    s0 * ( t0 * d0[coords(i0i, j0i)] +( t1 * d0[coords(i0i, j1i)]))
                   +s1 * ( t0 * d0[coords(i1i, j0i)] + ( t1 * d0[coords(i1i, j1i)]));
            }
        }
    set_bnd(b, d);
}
}
