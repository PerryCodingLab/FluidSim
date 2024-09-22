
final int SCALE = 5;
final int SPEED = 5;
FluidSqr fluid;
FluidSqr fluidTwo;
PGraphics canvas;
//int f = 1;



void settings() {
  size(N*SCALE, N*SCALE);
}

void setup(){
  //fluidTwo = new FluidSqr(1, 0.0005, 0.000005);
  fluid = new FluidSqr(2, 0.01, 0.0000005);
}

//void keyPressed(){
//  switch (key){
//    case '1': f = 1; break;
//    case '2': f = 2; break;
//    //case 'r': restart(); break;
//  }
//}

void mouseDragged(){
 fluid.addDensity(mouseX/SCALE, mouseY/SCALE, 5000);
 
 fluid.addVelocity(mouseX/SCALE, mouseY/SCALE, SPEED*(mouseX - pmouseX), SPEED*(mouseY-pmouseY));
}

void draw(){
  background(0);
  //stroke(51);
  //strokeWeight(2);
  //fluid.addDensity((N/2),(N/2),1000);
  //fluid.addVelocity((N/2),(N/2),5,-10);
  //if (f == 1) {
  //  fluid.step();
  //  fluid.renderD();
  //}
  //if (f == 2) {
  //  fluidTwo.step();
  //  fluidTwo.renderD();
  //}
  //fluidTwo.step();
  fluid.step();
  fluid.renderD();
  //fluidTwo.renderD();
}
