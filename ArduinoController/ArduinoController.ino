/* README 
  DO NOT USE MORE THAN 8 CHARACTERS TO SAVE THE FILE. SHIT!!!!!!!!!!! 
  READY TO FLYYYYYYYYYYYYYY


*/ 



#include <Servo.h>
#include <SD.h>
#include <Adafruit_APDS9960.h>
#include <Adafruit_BMP280.h>
#include <Adafruit_LIS3MDL.h>
#include <Adafruit_LSM6DS33.h>
#include <Adafruit_LSM6DS3TRC.h>
#include <Adafruit_SHT31.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_NeoPixel.h>    //  Library that provides NeoPixel functions

int Use_Simulated_Altitude  =0; 
///////////////////////////////////// TIME PARAMS /////////////////////////// 
int nextSensorUpdate; 
const int SensorUpdateInterval =200; // interval to read data from barometric pressure 
const int PIDUpdateInterval = 20000;
const int DataLoggerOutputInterval = 1000;
int nextDataLoggerOutput; 
int nextPIDUpdate;
int servoCloseTime;

/////////////////////////////////////////////////////////////////////////////// 

////////////////// THIS IS FOR DATALOGGING //////////////////////////// 
bool isFileInitialized = false; // To track if headers have been written
File myFile;
/////////////////////////////////////////////////////////////////////////////////////// 

////////////////// THIS IS FOR pressure sensor //////////////////////////// 
Adafruit_BMP280 bmp280;     // temperautre, barometric pressure 
/////////////////////////////////////////////////////////////////////////////////////// 


//////////////////////////// pre-defined the functions//////////////////////// 
float GetAltitudeFromSensor (); 
float GetAltitudeFromSimulation(String receivedString);
/////////////////////////////////////////////////////////////////////////////// 

/////////////////////////////// PARAMS FOR PID ////////////////////////////////
float CurrentAltitude;
// VARIABLES AND LIBRARY FOR PID 
float velocity_sp;  // The setpoint velocity
float avg_velocity =7;  // Average velocity
float estimated_altitude;  // Estimated altitude
float altitude_sp_min = 10000 ; // Minimum altitude setpoint
float altitude_sp = 19000;
float altitude_sp_max = 30000;
float time_constant = (altitude_sp-altitude_sp_min)/avg_velocity ;  // Time constant
float u =0;
float Kp = 0.4;
float Ki = 0;
float Kd = 0.1;
float error = 0;
float error_prev = 0;
float sum_int =0;
float control_interval = 20;
float d_error = 0;
String statusString;
float target_velocity =2;  // Target velocity
float controlParameter;
float pressure ;

/////////////////////////////////////////////////////////////////////////////////////// 


/////////////////// PARAMS FOR Linear Regression ///////////////////////// 
float S_xy =0.0;
float S_xx = 0.0;
float S_x = 0.0;
float S_y = 0.0; 
float Sum_xx =0.0; 
float Sum_xy = 0.0; 
float velocity = 0.0;
int datapoints =0;
float x,y;

//////////////////////// PARAMS FOR SERVO MOTOR ///////////// 
bool isValveOpen = false;
float ControlParameterMin = 0.15;
const int Degree_OPEN = 60;
const int Degree_CLOSE = 0;

Servo servo_9;

int currentTime;
Adafruit_NeoPixel onePixel = Adafruit_NeoPixel(1, 8, NEO_GRB + NEO_KHZ800);

void setup (){

  Serial.begin (115200);
  // initialize the sensors
  bmp280.begin();                                                                 

  // Initialize the build-in LED for debugging 
  onePixel.begin();             // Start the NeoPixel object
  onePixel.clear();             // Set NeoPixel color to black (0,0,0)
  onePixel.setBrightness(20);   // Affects all subsequent settings
  onePixel.show();              // Update the pixel state

  // Initialize SD card for datalogging 
  SD.begin (10);

  // Initialize servo - start from close position 
  servo_9.attach(9); // Attach servo to pin 9
  servo_9.write (0);
  delay (5000);
}


void loop() {

  currentTime = millis(); // use to decide which actions to take
  if (Use_Simulated_Altitude ==1){
      /////////////////////////////// READ DATA FROM SIMULATION /////////////////////////////////// UNCOMMENT ALL THE CODE INSIDE TO USE WITH SIMULATION 
      // if there is data in the buffer, read data from simulation
      while(Serial.available() >0){
        String receivedString = Serial.readStringUntil ('\n'); 
        CurrentAltitude = GetAltitudeFromSimulation (receivedString);
        // Calculate accumulator 
        x = currentTime/1000.0 ;
        y = CurrentAltitude;
        S_xy += x*y;
        S_x += x ;
        S_y +=y ;
        S_xx += x*x ;
        datapoints+=1;
      }
    }
  ////////////////////////////////////// Get Data  //////////////////////////// 
  if (currentTime > nextSensorUpdate) {
    /* 
      this is probably on a 200 ms interval
      read value from sensor and add to the accumulators for averaging and linear regression
    */ 
    if (Use_Simulated_Altitude ==0){
      //////////////////////////////////// READ DATA FROM SENSOR /////////////////////////////////// UNCOMMENT ALL THE CODE INSIDE TO USE WITH REAL SYSTEM
      
      CurrentAltitude = GetAltitudeFromSensor (); 
      pressure = bmp280.readPressure();

      if (true) { // this is to avoid when the sensor went wrong 
        // Calculate accumulator 
          x = currentTime/1000.0 ;
          y = CurrentAltitude;
          S_xy += x*y;
          S_x += x ;
          S_y +=y ;
          S_xx += x*x ;
          datapoints+=1;
      }
      else {
        // DO NOTTHING 
      }
    }
    
   // when you are using the simulator to generate this data, you will need to read each simulated altitude
   // data point from the serial stream – this section of code will be different when you are using the barometric
   // altimeter.
    //////////////////////////////////////////////////////////////////////////////////////////////
    // when using the simulator, you may need to be careful about timing of how often the simulator sends data,
    // and it is OK to miss data points (i.e. don’t always need 100 data points during every 20 second period)
    nextSensorUpdate = currentTime + SensorUpdateInterval;
  
  }

  if (currentTime > nextPIDUpdate) {

    ///////////////////////// CALCULATE AVERAGE ASCEND RATE DURING THE PRIOR CONTROL_INTERVAL///////////////////////
      if (datapoints!=0){ // in case we do not receive any data => prevent divide by zero 
        Sum_xx = S_xx - (S_x*S_x)/datapoints ; // change to N to see the big result 
        Sum_xy = S_xy - (S_x * S_y)/datapoints;
        velocity = Sum_xy/Sum_xx;
        estimated_altitude = CurrentAltitude; 
      }
      else { // IF WE DO NOT RECEIVE ANY DATAPOINTS => VELOCITY = VELOCITY_PREVIOUS 
        velocity = velocity;
        // what to do when sensor fails 
      }

    /////////////////////////////////////////////////////// PID ALGORITHM ////////////////////////////////////////// 
    // PID CONTROLLER 
      // GREEN MODE 
      if (estimated_altitude <altitude_sp_min){
        u =0;
        velocity_sp = 0; // infinity ? 
      } 
      // YELLOW MODE 
      else if ((altitude_sp_min <=estimated_altitude)  &&  (estimated_altitude < altitude_sp)){

        velocity_sp = avg_velocity - abs((estimated_altitude - altitude_sp_min) / time_constant) + target_velocity;
        // CALCULATE CONTROL SIGNAL 
        error = velocity-velocity_sp ;
        sum_int += error * control_interval ;
        d_error = (error - error_prev) / control_interval ;
        u =  Kp * error + Ki * sum_int + Kd * d_error;
        u = max(0.0, min(1.0, u));
        error_prev = error;
      }
      // RED MODE 
      else {
        velocity_sp =0;
        // CALCULATE CONTROL SIGNAL 
        error = velocity-velocity_sp ;
        sum_int += error * control_interval ;
        d_error = (error - error_prev) / control_interval ;
        u =  Kp * error + Ki * sum_int + Kd * d_error;
        u = max(0.0, min(1.0, u));
        error_prev = error;
      }
      // update the control signal to run the servo motor 
      controlParameter =u;
      // If use simulation => send control signal and velocity to Python
      if (Use_Simulated_Altitude==1){
        SendControlSignalToSimulation (u,velocity);
      }

    // DONE FOR CURRENT CONTROL LOOP => RESET ACCUMULATORS FOR NEXT CONTROL LOOP
      S_xy = S_x = S_y = S_xx = Sum_xx = Sum_xy = 0.0;
      datapoints =0;
    // OPEN SERVO MOTOR IF CONTROL PARAMETER IS VALID 
      if (controlParameter > ControlParameterMin){
        servo_9.write(Degree_OPEN); // open the servo motor 
        isValveOpen = true; // update the status of the servo motor 
        // delay (300); 
      }
    // INCREMENT TIME PARAMETERS 
      nextPIDUpdate = currentTime + PIDUpdateInterval;
      servoCloseTime = int (currentTime + float(PIDUpdateInterval * controlParameter)); 
        // NOTE control Parameter may be a float – so need to make this line work properly // I DO NOT UNDERSTAND WHAT DO YOU MEAN BY FLOAT 

    }

  if (currentTime > nextDataLoggerOutput) { 
    DataLogging (u,currentTime,CurrentAltitude,velocity,pressure); // control, time, altitude, velocity, pressure
    nextDataLoggerOutput = currentTime + DataLoggerOutputInterval;

    /* this is probably on a 1-2 second interval
      print out the current values of many variables:
      time, altitude, ascent rate, control parameter, vent state (open or closed), temperature, etc.
    */
  }
  /////////////////////////////////////// CLOSE THE SERVO MOTOR AT WAITING FOR TIME-INTERVAL //////////////
  if (currentTime > servoCloseTime) {
     // if servo is open, close servo -> set servo pin to close position value
     if (isValveOpen){
        servo_9.write(Degree_CLOSE); // close the servo motor 
        isValveOpen = false; // update the status of the servo motor 
     }
     
    // with our current vent control strategy, this section of code does not need to be complicated to handle timing
    // if the PID indicates that the servo should open during the control interval, then the PID should open the servo
    // then it only needs to be closed once – whenever (currentTime > nextServoCloseTime)
  }

} // end of loop function
 
float GetAltitudeFromSensor (){
  int altitude = bmp280.readAltitude(1013.25);
  return altitude ;
}
////////////////////////////////// GOOD TO GO ////////////////////////////////////////////////////
float GetAltitudeFromSimulation (String receivedString){
  float AltitudeFromSimulation;
  /* 
  This function return the altitude from the simulation
  This function is copied from the accelerated controller - SHOULD NOT HAVE ANYTHING PROBLEMS WITH RECEIVE DATA FROM PYTHON 
  */ 

  // Calculate the number of signals received 
    if (receivedString.length() > 0) {

  // Find the position of the semicolon
    int separatorIndex = receivedString.indexOf(';');
      if (separatorIndex > 0) { // Ensure the separator is found
        // Extract the altitude and time substrings
        String altitudeString = receivedString.substring(0, separatorIndex);
        // String timeString = receivedString.substring(separatorIndex + 1); // DO NOT CARE 
        // Convert the strings to floats
         AltitudeFromSimulation = altitudeString.toFloat();
        // time_i = timeString.toFloat();  // DO NOT CARE ABOUT THIS  
 
        // delay (time_delay);
      }
    }
  
  return AltitudeFromSimulation;
} 

void SendControlSignalToSimulation (float u, float velocity){
  Serial.print (u,2);
  Serial.print(";");
  Serial.print (velocity,2);
  Serial.println();
  // delay  (time_delay);
}


// THis function will be put into the for loop 
void DataLogging (float u_save, float time_Arduino_save,float altitude_save, float velocity_save,float pressure_save){
  /* 
  This function receives the data, then save to SD card 
  DO NOT USE MORE THAN 8 CHARACTERS TO SAVE THE FILE. SHIT!!!!!!!!!!! 

  */ 
  
  // Open the file in append mode
  myFile = SD.open("test4.csv", FILE_WRITE);  // DO NOT USE MORE THAN 8 CHARACTERS TO SAVE THE FILE. SHIT!!!!!!!!!!! 

  
  if (!myFile) {
    Serial.println("Error opening file for writing!");
    return;
  }

  if (myFile) {
    // Write headers if the file is empty or just created
    if (!isFileInitialized) {
      myFile.println("Time Arduino (s), Altitude (m), Control signal, Velocity (m/s), Pressure (Pa)");
      isFileInitialized = true; // Mark headers as written
    }

    // Write data as CSV format
    myFile.print(time_Arduino_save, 2);
    myFile.print(",");
    myFile.print(altitude_save);
    myFile.print(",");
    myFile.print(u_save);
    myFile.print(",");
    myFile.print(velocity_save);
    myFile.print(",");
    myFile.println(pressure_save);
    myFile.close(); // Close the file after writing
    onePixel.setPixelColor(0, 0, 255, 0);   // GREEN for success
  } 
  else {
    onePixel.setPixelColor(0, 255, 0, 0);   // RED for failure
  }

  onePixel.show();      
  delay(50);
  onePixel.clear(); 
  onePixel.show();
  delay (50);
}





