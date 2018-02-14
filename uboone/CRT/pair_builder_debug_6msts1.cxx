///////////////////////////////////////////////////////////////////////////////////////////////////////
//    This program receive hits from a febdriver and searches for ts1_ref events
//    Note: Set maximal timedifference of events to select wrt ts1_ref in the code
//    This Version is for permanent use (no counters)
//    This filter is for uBooNe!!!
//    Written by Thomas Mettler
///////////////////////////////////////////////////////////////////////////////////////////////////////


#include "pair_builder_debug_6msts1.h"
using namespace crt;

crt::EVENT_t evbuf[MAXFEBNR*EVSPERFEB+1];    //buffer to receive events (same structure as the receiving events)
crt::EVENT_tpro evbuf_pro[MAXFEBNR+1][4*EVSPERFEB+1];  //buffer for processing (add the second, millisecond from sepcial events)
crt::EVENT_tpro evbuf_scan[MAXFEBNR+1][4*EVSPERFEB+1]; //buffer for scanning for coincidences (same structure as the buffer for processing)
crt::EVENT_tpro evbuf_filter[MAXFEBNR+1][4*EVSPERFEB+1];
crt::EVENT_tpro evbuf_filter_scan[MAXFEBNR+1][4*EVSPERFEB+1];
crt::EVENT_t_send beam_ev[10][4*EVSPERFEB+1];    //buffer to send out the coincidences (structure idealy same as the received events)
crt::EVENT_t ts0_ref_event[2];
crt::EVENT_t_send ts0_ref_event_buffer[MAXFEBNR+1][2];
crt::EVENT_t_send coincidence[10][MAXFEBNR+1];    //buffer to send out the coincidences (structure idealy same as the received events)
crt::EOP_EVENT_t refevent;
FILE *file_store;

std::vector<crt::CRTHit>  allmyCRTHits;

const static double inch=2.54; //inch in cm

void crt::pair_builder::make_pairs(const char * filename,const char * filename_store, int run_mode_, std::vector<crt::CRTHit>& allCRTHitstmp){
  InitReconst();
  find_pairs(filename, filename_store,run_mode_);
  allCRTHitstmp = allmyCRTHits;
}


int crt::pair_builder::find_pairs(const char * filename,const char * filename_store, int run_mode_){
 //if(argc!=4) { usage(); return 0;} //test the input variables and give hint if needed

  //char *filename=argv[1];
  //allmyCRTHits = allCRTHits;
  data=fopen(filename,"r");
  fseek(data, 0, SEEK_END); // seek to end of file
	long size = ftell(data); // get current file pointer
	fseek(data, 0, SEEK_SET); // seek back to beginning of file
	
	size_ev=size/sizeof(EVENT_t);		//number of total events
  //size_ev=100000;
	printf("Total Number of events: %ld\n",size_ev);

  //char *filename_store=argv[2];
  file_store=fopen(filename_store,"wb");
	run_mode=run_mode_;

  //set all hit counters (from scan buffer and pro buffer) to 0
  for(int i=0;i<MAXFEBNR+1;i++){
    ev_counter_mac[i]=0;
    ev_counter_scan[i]=0;
    ev_counter_filter[i]=0;
    order_buffer[i].flags=1;
    nrwm1[i]=0;
    nrwm2[i]=0;
    ts0ref_counter[i]=0;
    ev_counter_filter_scan[i]=0;
    event_time_diff[i]=0;
    event_time_diff_old[i]=0;
    previous_sec_scan[i]=0;
    previous_sec_mac[i]=0;
  }
  previous_sec=0;
  previous_ms=0;
  for(int i=0; i<100;i++){
    ts1ref_counter[i]=0;
    ts1ref_second[i]=0;
    send_bufnr=0;
    counter_tot=0;
    counter_old=0;
    stuck_event_counter=0;
    minus2_counter1=0;
    same_counter1=0;
    verbose_=0;
  }
  ready_to_fill=PROBUF_READY_TO_FILL;
  
  //endless loop for receiving-processing-sending events///////////////////////////////////////////////
  while(1){   //endless loop over all events receiving   
    //If one pro buffer is full->scan the whole buffer without scaling, else print status of buffer
    for(int i=0;i<MAXFEBNR;i++){
      if(ev_counter_mac[i]>(4*EVSPERFEB)){  //test if there is an overflow in the receiving buffer
        //error++; 
        printf("pro buffer scaned and reseted without scaling of %d...\n",i);
        shift_scale(i, 4*EVSPERFEB, 1e9);
        ready_to_fill=PROBUF_READY_TO_FILL;
        ev_counter_mac[i]=0;  //if one buffer is overload, it is scaled w/out scaling...
      }
      else if(ev_counter_mac[i]!=0){ //if everything if fine, print the number of events in the buffers
        //printf("fill status of %d: %d - %d\n",i,ev_counter_mac[i], ev_counter_scan[i]);
      }
    }
    //receive new data
    if(ready_to_fill==PROBUF_READY_TO_FILL){
      ready_to_fill=PROBUF_FILLING;
      if(counter_tot>size_ev) return 0;
      receive_data();}
    //scale and scan the new data
    if(ready_to_fill==PROBUF_READY_TO_SCALE){
      ready_to_fill=PROBUF_SHIFT_SCALE;
      scale_buffer();   //in case more than one second is in the buffer
      scale_buffer();   //scale_buffer should then be used more than ones
      //scale_buffer();   //only one is needed for polltimes<1 sec
    }
    else{ printf("pro buffer is in use... \n");
         //error++;
        }
    while(ready_to_fill!=PROBUF_READY_TO_SCALE && ready_to_fill!=PROBUF_READY_TO_FILL){
      //wait...
    }
  }
  fclose(file_store);
  return 0;
}

void crt::pair_builder::usage(){      //print some hints if input is wrong
 printf("Connects to data stream from running febdrv at a given data socket and send (adds) data.\n Usage: ");
 printf("beam_filter_file <path to file> <path_to_store> run_mode\n");
 printf("example:  \"rawdata_top.bin  file_store.bin 1\n");
}
//receive a poll of hits from the modules//////////////////////////////////////////////////////
void crt::pair_builder::receive_data(){
  int mac=0;
  for(int i=0;i<1000;i++){
    fread(&evbuf[i],sizeof(EVENT_t),1,data);
  }
  int received_events=1000;
  counter_tot+=received_events;
  if(counter_tot-counter_old>10000){
   printf("\rprocessed: %ld/%ld (%ld%%), stuckt events: %ld (%.2f permill), second assignement: same: %d, lost: %d", counter_tot, size_ev, 100*counter_tot/size_ev, stuck_event_counter, 1000*(float)stuck_event_counter/(counter_tot-1000), same_counter1,minus2_counter1);
    counter_old=counter_tot;
  }
  //int prereassign=0;//control number for changes in second assignement
  int reassign=0;//control number for changes in second assignement
  int check=0;//control number for print outs in second assignement
  int all_fine=0;//control number if the jump happened after a reference pulse
	
  for(int i=0; i<received_events; i++){
    //prereassign=0;
    if(evbuf[i].mac5==0xFFFF){    //reads the spezial event
      memcpy(&refevent,&evbuf[i].mac5,sizeof(EOP_EVENT_t));
      previous_sec=act_time[0][MAXFEBNR];
      previous_ms=act_time[1][MAXFEBNR];
      //printf("End:	%d seconds	%d millisec\n",act_time[0][MAXFEBNR], act_time[1][MAXFEBNR]);
      act_time[0][MAXFEBNR]=(int)refevent.end_s;
      act_time[1][MAXFEBNR]=(int)refevent.end_ms;
      //prozessed_hit_counter++;
      if(verbose_!=0){
        printf("%d\n",i);
      printf("pre+2: %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
      printf("start: sec: %10d, ts0: %10d, mac: %d, flags: %2d\n",act_time[0][mac],evbuf[i-1].ts0,evbuf[i-1].mac5, evbuf[i-1].flags);
      printf("start: sec: %10d, ts0: %10d, mac: %d, flags: %2d\n", evbuf_pro[evbuf[i-1].mac5][ev_counter_mac[evbuf[i-1].mac5]-1].sec,evbuf_pro[evbuf[i-1].mac5][ev_counter_mac[evbuf[i-1].mac5]-1].ts0,evbuf_pro[evbuf[i-1].mac5][ev_counter_mac[evbuf[i-1].mac5]-1].mac5, evbuf_pro[evbuf[i-1].mac5][ev_counter_mac[evbuf[i-1].mac5]-1].flags); 
      printf("Start: sec: %10d, mil: %10d millisec\n",(int)refevent.start_s, refevent.start_ms);
      printf("\n");
      printf("End:   sec: %10d, mil: %10d millisec\n",(int)refevent.end_s, refevent.end_ms);
      printf("\n");
      }
      for(int febnr=0;febnr<MAXFEBNR;febnr++){
        if(ev_counter_mac[febnr]>0){
          if(evbuf_pro[febnr][ev_counter_mac[febnr]-1].sec!=act_time[0][MAXFEBNR]){
            if((evbuf_pro[febnr][ev_counter_mac[febnr]-1].ts0<8e8) && (act_time[1][MAXFEBNR]<800) && (act_time[1][MAXFEBNR]>100)){
              //we have a problem here...
              //prereassign=1;
              if(verbose_!=0){
                printf("\nERROR!!!!      act-1=ev_sec    !!!!!!, mac= %d, lostcpu: %d\n", febnr, evbuf_pro[febnr][ev_counter_mac[febnr]-1].lostcpu );
                printf("start: sec: %10d, ts0: %10d, mac: %d, flags: %2d\n", evbuf_pro[febnr][ev_counter_mac[febnr]-1].sec,evbuf_pro[febnr][ev_counter_mac[febnr]-1].ts0,evbuf_pro[febnr][ev_counter_mac[febnr]-1].mac5, evbuf_pro[febnr][ev_counter_mac[febnr]-1].flags); 
                printf("End:   sec: %10d, mil: %10d millisec\n",(int)refevent.end_s, refevent.end_ms);
              }
              for(int corr=0; corr<ev_counter_mac[febnr]; corr++){
                evbuf_pro[febnr][corr].recover=11;
                evbuf_pro[febnr][corr].sec=act_time[0][MAXFEBNR];
                act_time[0][febnr]=act_time[0][MAXFEBNR];
              }
	      
            }
            else if((evbuf_pro[febnr][ev_counter_mac[febnr]-1].ts0>8e8) && (act_time[1][MAXFEBNR]<800) && (evbuf_pro[febnr][ev_counter_mac[febnr]-1].sec>act_time[0][MAXFEBNR])){
              //we have a big prolem here act+2=ev_sec
              if(verbose_!=0){
                printf("\nERROR!!!!      act+2=ev_sec    !!!!!!, mac= %d, lostcpu: %d\n", febnr, evbuf_pro[febnr][ev_counter_mac[febnr]-1].lostcpu );
                printf("start: sec: %10d, ts0: %10d, mac: %d, flags: %2d\n", evbuf_pro[febnr][ev_counter_mac[febnr]-1].sec,evbuf_pro[febnr][ev_counter_mac[febnr]-1].ts0,evbuf_pro[febnr][ev_counter_mac[febnr]-1].mac5, evbuf_pro[febnr][ev_counter_mac[febnr]-1].flags); 
                printf("End:   sec: %10d, mil: %10d millisec\n",(int)refevent.end_s, refevent.end_ms);
              }
              for(int corr=0; corr<ev_counter_mac[febnr]; corr++){
                //printf("start: sec: %10d, ts0: %10d, ts1: %10d, mac: %d, flags: %2d\n", evbuf_pro[febnr][corr].sec,evbuf_pro[febnr][corr].ts0,evbuf_pro[febnr][corr].ts1,evbuf_pro[febnr][corr].mac5, evbuf_pro[febnr][corr].flags); 
                evbuf_pro[febnr][corr].recover=12;
                evbuf_pro[febnr][corr].sec=act_time[0][MAXFEBNR];
                act_time[0][febnr]=act_time[0][MAXFEBNR];
              }
            }
            else if((evbuf_pro[febnr][ev_counter_mac[febnr]-1].ts0<8e8) && (act_time[1][MAXFEBNR]>800) && (evbuf_pro[febnr][ev_counter_mac[febnr]-1].sec<act_time[0][MAXFEBNR])){
              //we have a big prolem here act-2=ev_sec
              if(verbose_!=0){
                printf("\nERROR!!!!      act-2=ev_sec    !!!!!!, mac= %d, lostcpu: %d\n", febnr, evbuf_pro[febnr][ev_counter_mac[febnr]-1].lostcpu );
                printf("start: sec: %10d, ts0: %10d, mac: %d, flags: %2d\n", evbuf_pro[febnr][ev_counter_mac[febnr]-1].sec,evbuf_pro[febnr][ev_counter_mac[febnr]-1].ts0,evbuf_pro[febnr][ev_counter_mac[febnr]-1].mac5, evbuf_pro[febnr][ev_counter_mac[febnr]-1].flags); 
                printf("End:   sec: %10d, mil: %10d millisec\n",(int)refevent.end_s, refevent.end_ms);
              }
              for(int corr=0; corr<ev_counter_mac[febnr]; corr++){
                //printf("start: sec: %10d, ts0: %10d, ts1: %10d, mac: %d, flags: %2d\n", evbuf_pro[febnr][corr].sec,evbuf_pro[febnr][corr].ts0,evbuf_pro[febnr][corr].ts1,evbuf_pro[febnr][corr].mac5, evbuf_pro[febnr][corr].flags); 
                evbuf_pro[febnr][corr].recover=13;
                evbuf_pro[febnr][corr].sec=act_time[0][MAXFEBNR];
                act_time[0][febnr]=act_time[0][MAXFEBNR];
              }
            }
          } 
        } 
      }
    }
    else {
      //add the second to the events and copy the buffer into another
      mac=evbuf[i].mac5; 
      check=0;
      reassign=0;
      if(ev_counter_mac[mac]>0 && (evbuf_pro[mac][ev_counter_mac[mac]-1].flags==3 || evbuf_pro[mac][ev_counter_mac[mac]-1].flags==1 || evbuf_pro[mac][ev_counter_mac[mac]-1].lostcpu==99)){ //check for jump, excluding all ref_events
        if(my_abs(event_time_diff[mac],(evbuf[i].ts0-evbuf[i].ts1))>500){ // the jump ha to be bigger than 500
          all_fine=0;
          //check for what referent event was missed or if it was a stucked event...
          if(my_abs(event_time_diff[mac],(evbuf[i].ts0-evbuf[i].ts1)-1e9)<30000){ //check if ts0_ref was missed 30us > dead time FEB
            evbuf[i].lostcpu=1;
            if(evbuf_pro[mac][ev_counter_mac[mac]-1].flags==5 || evbuf_pro[mac][ev_counter_mac[mac]-1].flags==7) all_fine=1;
            if(all_fine!=1){
              if(previous_sec==act_time[0][MAXFEBNR]) act_time[0][mac]=act_time[0][MAXFEBNR]+1; //check if already poll of new second
              else act_time[0][mac]=act_time[0][MAXFEBNR];
              if(act_time[0][mac]==previous_sec_mac[mac]){reassign=1; check=2;}
              if(act_time[0][mac]==previous_sec_mac[mac]+2){
                //printf("1pre+2 %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
                //act_time[0][mac]--;
                reassign=3;
                check=2;
              }
              if(evbuf_pro[mac][ev_counter_mac[mac]-1].sec==act_time[0][mac] && reassign!=1) {same_counter1++;check=1;}
              if(my_abs(evbuf_pro[mac][ev_counter_mac[mac]-1].sec,act_time[0][mac])>1&& reassign!=3) {
                //printf("1pre+2.2 %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
                reassign=3;
                act_time[0][mac]--;
                minus2_counter1++;
                check=2;}
              ts0ref_counter[mac]++;
              //assign new event at the right place with ts0 = 1e9 (no scaling)
              evbuf_pro[mac][ev_counter_mac[mac]].sec=act_time[0][mac];
              evbuf_pro[mac][ev_counter_mac[mac]].ms=act_time[1][MAXFEBNR];
              evbuf_pro[mac][ev_counter_mac[mac]].mac5=mac;
              evbuf_pro[mac][ev_counter_mac[mac]].flags=7;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0=1e9;
              evbuf_pro[mac][ev_counter_mac[mac]].ts1=evbuf_pro[mac][ev_counter_mac[mac]-1].ts1+(1e9-evbuf_pro[mac][ev_counter_mac[mac]-1].ts0);
              for(int j=0; j<32;j++) evbuf_pro[mac][ev_counter_mac[mac]].adc[j]=0;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0_scaled=0;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0_ref=0;
              evbuf_pro[mac][ev_counter_mac[mac]].lostcpu=0;
              evbuf_pro[mac][ev_counter_mac[mac]].lostfpga=0;
              //printf("Found TS0_ref!!! new event: mac: %d, flags: %2d, ts0: %10d, ts1: %10d, sec: %10d\n", evbuf_pro[mac][ev_counter_mac[mac]].mac5, evbuf_pro[mac][ev_counter_mac[mac]].flags, evbuf_pro[mac][ev_counter_mac[mac]].ts0,evbuf_pro[mac][ev_counter_mac[mac]].ts1,evbuf_pro[mac][ev_counter_mac[mac]].sec); 
              ev_counter_mac[mac]++;
            }
          }
          else if(my_abs(event_time_diff_old[mac],(evbuf[i].ts0+1e9-evbuf[i].ts1))<30000){ //check if ts0_ref was missed 30us > dead time FEB
            evbuf[i].lostcpu=2;
            if(evbuf_pro[mac][ev_counter_mac[mac]-1].flags==5 || evbuf_pro[mac][ev_counter_mac[mac]-1].flags==7) all_fine=1;
            if(all_fine!=1){
              if(previous_sec==act_time[0][MAXFEBNR]) act_time[0][mac]=act_time[0][MAXFEBNR]+1; //check if already poll of new second
              else act_time[0][mac]=act_time[0][MAXFEBNR];
              if(act_time[0][mac]==previous_sec_mac[mac]){reassign=1; check=2;}
              if(act_time[0][mac]==previous_sec_mac[mac]+2){
                //printf("2pre+2 %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
                //act_time[0][mac]--;
                reassign=3;
                check=2;
              }
              if(evbuf_pro[mac][ev_counter_mac[mac]-1].sec==act_time[0][mac] && reassign!=1) {same_counter1++;check=1;}
              if(my_abs(evbuf_pro[mac][ev_counter_mac[mac]-1].sec,act_time[0][mac])>1&& reassign!=3) {
                //printf("2pre+2.2 %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
                reassign=3;
                act_time[0][mac]--;
                minus2_counter1++;
                check=2;}
              ts0ref_counter[mac]++;
              //assign new event at the right place with ts0 = 1e9 (no scaling)
              evbuf_pro[mac][ev_counter_mac[mac]].sec=act_time[0][mac];
              evbuf_pro[mac][ev_counter_mac[mac]].ms=act_time[1][MAXFEBNR];
              evbuf_pro[mac][ev_counter_mac[mac]].mac5=mac;
              evbuf_pro[mac][ev_counter_mac[mac]].flags=7;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0=1e9;
              evbuf_pro[mac][ev_counter_mac[mac]].ts1=evbuf_pro[mac][ev_counter_mac[mac]-1].ts1+(1e9-evbuf_pro[mac][ev_counter_mac[mac]-1].ts0);
              for(int j=0; j<32;j++) evbuf_pro[mac][ev_counter_mac[mac]].adc[j]=0;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0_scaled=0;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0_ref=0;
              evbuf_pro[mac][ev_counter_mac[mac]].lostcpu=0;
              evbuf_pro[mac][ev_counter_mac[mac]].lostfpga=0;
              //printf("Found TS0_ref!!! new event: mac: %d, flags: %2d, ts0: %10d, ts1: %10d, sec: %10d\n", evbuf_pro[mac][ev_counter_mac[mac]].mac5, evbuf_pro[mac][ev_counter_mac[mac]].flags, evbuf_pro[mac][ev_counter_mac[mac]].ts0,evbuf_pro[mac][ev_counter_mac[mac]].ts1,evbuf_pro[mac][ev_counter_mac[mac]].sec); 
              ev_counter_mac[mac]++;
            }
          }
          else if((evbuf[i].ts1<(evbuf[i].ts0-evbuf_pro[mac][ev_counter_mac[mac]-1].ts0+20)&& evbuf[i].ts0<(evbuf_pro[mac][ev_counter_mac[mac]-1].ts0+1e8))){ //check if ts1_ref 
            evbuf[i].lostcpu=3;
            if(evbuf_pro[mac][ev_counter_mac[mac]-1].flags==10 || evbuf_pro[mac][ev_counter_mac[mac]-1].flags==11) all_fine=1;
            if(all_fine!=1){
              if(my_abs(evbuf[i].ts0,(evbuf_pro[mac][ev_counter_mac[mac]-1].ts0-evbuf[i].ts1))<20) evbuf_pro[mac][ev_counter_mac[mac]-1].flags=11; //1. check if there is an event with wrong assigned flag
              else if((evbuf[i].ts0-evbuf_pro[mac][ev_counter_mac[mac]-1].ts0)<1e7){  //2. if no event is there add new event
              evbuf_pro[mac][ev_counter_mac[mac]].sec=act_time[0][mac];
              evbuf_pro[mac][ev_counter_mac[mac]].ms=act_time[1][MAXFEBNR];
              evbuf_pro[mac][ev_counter_mac[mac]].mac5=mac;
              evbuf_pro[mac][ev_counter_mac[mac]].flags=11;//evbuf[i].flags | 0x100;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0=evbuf[i].ts0-evbuf[i].ts1;
              evbuf_pro[mac][ev_counter_mac[mac]].ts1=evbuf_pro[mac][ev_counter_mac[mac]-1].ts1+(evbuf[i].ts0-evbuf[i].ts1-evbuf_pro[mac][ev_counter_mac[mac]-1].ts0);
              for(int j=0; j<32;j++) evbuf_pro[mac][ev_counter_mac[mac]].adc[j]=0;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0_scaled=0;
              evbuf_pro[mac][ev_counter_mac[mac]].ts0_ref=0;    //not really used
              evbuf_pro[mac][ev_counter_mac[mac]].lostcpu=0;
              evbuf_pro[mac][ev_counter_mac[mac]].lostfpga=0;
              //printf("\nnew event: mac: %d, flags: %2d, ts0: %10d, ts1: %10d, sec: %10d\n", evbuf_pro[mac][ev_counter_mac[mac]].mac5, evbuf_pro[mac][ev_counter_mac[mac]].flags, evbuf_pro[mac][ev_counter_mac[mac]].ts0,evbuf_pro[mac][ev_counter_mac[mac]].ts1,evbuf_pro[mac][ev_counter_mac[mac]].sec); 
              ev_counter_mac[mac]++;}
              else {evbuf[i].lostcpu=99; stuck_event_counter++;}
            }
          }
          else {  //it has to be a stucked event. Count them for statistic, asign something special to it...
            evbuf[i].lostcpu=99;
            //mac=1;
            stuck_event_counter++;
          }
        }
        else{
          if(evbuf_pro[mac][ev_counter_mac[mac]-1].lostcpu==99) evbuf[i].lostcpu=10;
        }
        //event_time_diff[mac]=evbuf[i].ts0-evbuf[i].ts1; //calculating again the timedifference between ts0 and ts1, for jump searches in the beginning...
      }

      //fill the aktual event normal (if a jump happend, an event before is inserted if a ref was missed...)
      if((evbuf[i].flags==7 || evbuf[i].flags==5) && (evbuf[i].lostcpu!=99&&evbuf[i].lostcpu!=10)){
        //printf("TS=refeventstart: sec: %10d, ts0: %10d, mac: %d, flags: %2d\n",act_time[0][mac],evbuf[i].ts0,evbuf[i].mac5, evbuf[i].flags);
        previous_sec_mac2[mac]=previous_sec_mac[mac];
        previous_sec_mac[mac]=act_time[0][mac];
        if(previous_sec==act_time[0][MAXFEBNR]) act_time[0][mac]=act_time[0][MAXFEBNR]+1;
        else act_time[0][mac]=act_time[0][MAXFEBNR];
        if(act_time[0][mac]==previous_sec_mac[mac]){reassign=1; check=2;}
        if(act_time[0][mac]==previous_sec_mac[mac]+2){
          //printf("pre+2: %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
          if(act_time[1][MAXFEBNR]>200) act_time[0][mac]--;
          reassign=2;
          check=2;
        }
        //act_time[0][mac]=act_time[0][MAXFEBNR]+1;  // || evbuf[i].ts0<=(act_time[1][mac])
        ts0ref_counter[mac]++;
        if(evbuf_pro[mac][ev_counter_mac[mac]-1].sec==act_time[0][mac] && reassign!=1){
          same_counter1++;
          //printf("same %d, %d, %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac],evbuf_pro[mac][ev_counter_mac[mac]-1].sec,evbuf_pro[mac][ev_counter_mac[mac]-1].ts0);
          check=1;
        }
        if(my_abs(evbuf_pro[mac][ev_counter_mac[mac]-1].sec,act_time[0][mac])>1&& reassign!=2) {
          //printf("pre+2.2 %d, %d, %d\n",previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
          reassign=2;
          act_time[0][mac]--;
          minus2_counter1++;
          check=2;}
        //act_time[1][mac]=act_time[1][MAXFEBNR];
        //printf("ts0: %d, ts1: %d\n", first, second);
      }
	    if((evbuf[i].flags==7 || evbuf[i].flags==5)&&(evbuf[i].lostcpu==99||evbuf[i].lostcpu==10)) ts0ref_counter[mac]++;

      evbuf_pro[mac][ev_counter_mac[mac]].sec=act_time[0][mac];
      evbuf_pro[mac][ev_counter_mac[mac]].ms=act_time[1][MAXFEBNR];
      evbuf_pro[mac][ev_counter_mac[mac]].mac5=mac;
      evbuf_pro[mac][ev_counter_mac[mac]].flags=evbuf[i].flags;
      evbuf_pro[mac][ev_counter_mac[mac]].ts0=evbuf[i].ts0;
      evbuf_pro[mac][ev_counter_mac[mac]].ts1=evbuf[i].ts1;
      for(int j=0; j<32;j++) evbuf_pro[mac][ev_counter_mac[mac]].adc[j]=evbuf[i].adc[j];
      evbuf_pro[mac][ev_counter_mac[mac]].ts0_scaled=0;
      evbuf_pro[mac][ev_counter_mac[mac]].ts0_ref=0;    //not really used
	    evbuf_pro[mac][ev_counter_mac[mac]].recover=0;  
      evbuf_pro[mac][ev_counter_mac[mac]].lostcpu=evbuf[i].lostcpu;
      evbuf_pro[mac][ev_counter_mac[mac]].lostfpga=evbuf[i].lostfpga;
      evbuf_pro[mac][ev_counter_mac[mac]].nrtrigger_11=ts0ref_counter[mac];
      //printf("End:   sec  %d ms: %d\n",(int)refevent.end_s, refevent.end_ms);
      //printf("event: sec: %d ms: %d ....%d mac5,  ts0: %d \n", evbuf_pro[mac][ev_counter_mac[mac]].sec, evbuf_pro[mac][ev_counter_mac[mac]].ms,evbuf_pro[mac][ev_counter_mac[mac]].mac5, evbuf_pro[mac][ev_counter_mac[mac]].ts0);
      ev_counter_mac[mac]++;
      if(reassign==1){
        for(int z=0; z<ev_counter_mac[mac]-1;z++){
          evbuf_pro[mac][z].sec=evbuf_pro[mac][z].sec-1;
		      evbuf_pro[mac][z].recover=1;
        }
        check=1;
      }
      if(reassign==2){
        if(act_time[1][MAXFEBNR]<200){
            for(int z=0; z<ev_counter_mac[mac]-1;z++){
              evbuf_pro[mac][z].sec=evbuf_pro[mac][z].sec+1;
              evbuf_pro[mac][z].recover=2;
            }
            if(evbuf_pro[mac][ev_counter_mac[mac]].flags!=5 || evbuf_pro[mac][ev_counter_mac[mac]].flags!=7) evbuf_pro[mac][ev_counter_mac[mac]].sec+=1;
        }
        //else act_time[0][mac]--;
		    check=2;
      }
      
      if(reassign==3){
        for(int z=0; z<ev_counter_mac[mac]-1;z++){
         //evbuf_pro[mac][z].sec=evbuf_pro[mac][z].sec+1;
          evbuf_pro[mac][z].recover=3;
        }
        //if(evbuf_pro[mac][ev_counter_mac[mac]].flags!=5 || evbuf_pro[mac][ev_counter_mac[mac]].flags!=7) evbuf_pro[mac][ev_counter_mac[mac]].sec+=1;
      }
      
      if(verbose_!=0){ //just printouts if needed... Row of events read in...
        printf("check: %d\n", check);
        printf("End:	%d seconds	%d millisec\n",act_time[0][MAXFEBNR], act_time[1][MAXFEBNR]);
		    printf("pre: %d: %d, %d, %d\n",reassign,previous_sec_mac2[mac], previous_sec_mac[mac], act_time[0][mac]);
        int sum=0, sum1=0, sum2=0, sum3=0, sum4=0;
        for(int amp=0; amp<32;amp++){
          sum+=evbuf[i].adc[amp];
          sum1+=evbuf_pro[mac][ev_counter_mac[mac]].adc[amp];
          sum3+=evbuf_pro[mac][ev_counter_mac[mac]-1].adc[amp];
          sum4+=evbuf_pro[mac][ev_counter_mac[mac]-2].adc[amp];
          sum2+=evbuf[i+1].adc[amp];
        }
          //printf("after:  previous sec: %d, this second: %d\n",previous_sec,act_time[0][MAXFEBNR]);
        for(int z=0; z<5;z++){
          //printf("jump: %3d: mac: %2d, flags: %2d, ts0: %10d, ts1: %10d, sec: %10d  dt: %d, lostcpu: %d\n", z,evbuf_pro[mac][z].mac5, evbuf_pro[mac][z].flags, evbuf_pro[mac][z].ts0,evbuf_pro[mac][z].ts1,evbuf_pro[mac][z].sec, evbuf_pro[mac][z].ts0-evbuf_pro[mac][z].ts1, evbuf_pro[mac][z].lostcpu ); 
          //evbuf_pro[mac][z].sec--;
         }
        for(int z=ev_counter_mac[mac]-5; z<ev_counter_mac[mac];z++){
          //for(int z=0; z<ev_counter_mac[mac];z++){
          printf("jump: %3d: mac: %2d, flags: %2d, ts0: %10d, ts1: %10d, sec: %10d  dt: %d, lostcpu: %d\n", z,evbuf_pro[mac][z].mac5, evbuf_pro[mac][z].flags, evbuf_pro[mac][z].ts0,evbuf_pro[mac][z].ts1,evbuf_pro[mac][z].sec, evbuf_pro[mac][z].ts0-evbuf_pro[mac][z].ts1, evbuf_pro[mac][z].lostcpu ); 
          //evbuf_pro[mac][z].sec--;
         }
        printf("jump+1: mac: %d, flags: %2d, ts0: %10d, ts1: %10d, sum adc: %5d sec: %10d, dt: %d\n", evbuf[i+1].mac5, evbuf[i+1].flags, evbuf[i+1].ts0,evbuf[i+1].ts1, sum2,0, evbuf[i+1].ts0-evbuf[i+1].ts1); 
        printf("jump+2: mac: %d, flags: %2d, ts0: %10d, ts1: %10d, sum adc: %5d sec: %10d, dt: %d\n", evbuf[i+2].mac5, evbuf[i+2].flags, evbuf[i+2].ts0,evbuf[i+2].ts1, sum2,0, evbuf[i+2].ts0-evbuf[i+2].ts1); 
        printf("jump+3: mac: %d, flags: %2d, ts0: %10d, ts1: %10d, sum adc: %5d sec: %10d, dt: %d\n", evbuf[i+3].mac5, evbuf[i+3].flags, evbuf[i+3].ts0,evbuf[i+3].ts1, sum2,0,evbuf[i+3].ts0-evbuf[i+3].ts1); 
        printf("jump+4: mac: %d, flags: %2d, ts0: %10d, ts1: %10d, sum adc: %5d sec: %10d, dt: %d\n", evbuf[i+4].mac5, evbuf[i+4].flags, evbuf[i+4].ts0,evbuf[i+4].ts1, sum2,0,evbuf[i+4].ts0-evbuf[i+4].ts1);
        printf("\n");
      }
      if(evbuf[i].flags!=99) { //if the last one was not a stucked event, calculate the ts0-ts1 new...
        if(abs(event_time_diff[mac]-event_time_diff_old[mac])>50) event_time_diff_old[mac]=event_time_diff[mac];
        event_time_diff[mac]=evbuf[i].ts0-evbuf[i].ts1;
      }
    }
  }
  ready_to_fill=PROBUF_READY_TO_SCALE;
}
//searches referent events (produced throw PPS) and process them further/////////////////////////////
void crt::pair_builder::scale_buffer(){
  //printf("scale all buffers\n");
  for(int j=0;j<MAXFEBNR;j++){ //loop over all planes
   for(int i=0;i<ev_counter_mac[j];i++){  
    if(evbuf_pro[j][i].flags==5 || evbuf_pro[j][i].flags==7){ //search a time ref event
      if(ready_to_scan==SCANBUF_READY_TO_FILL){
        order_buffer[j].sec=evbuf_pro[j][i].sec;
        order_buffer[j].ref_nr=i;
        order_buffer[j].ts0_ref=evbuf_pro[j][i].ts0;
        order_buffer[j].flags=0;
        //printf("%d %d %d %d %d %d\n", evbuf_pro[j][i].mac5, evbuf_pro[j][i].flags,evbuf_pro[j][i].ts0, evbuf_pro[j][i].ts1, evbuf_pro[j][i].adc[0],evbuf_pro[j][i].sec);
      }
      else {printf("Buffer not ready to scan!!\n"); 
            //error++;
           }
      //printf("break!\n");
      i=ev_counter_mac[j];  //to go of the for loop
      break;
    }
   }   
  }
  //the following part looks that the buffers with sec<sec_max is scaned first/////
  uint32_t sec_max=0;
  for(int j=0;j<MAXFEBNR+1;j++){  //find max_sec
    if(order_buffer[j].flags==0 && order_buffer[j].sec>sec_max) sec_max=order_buffer[j].sec;
  } 
  //Scan the buffer first with second < max_second in the buffer
  for(int j=0;j<MAXFEBNR+1;j++){  //scale/scan all seconds < max_sec
    if(order_buffer[j].flags==0 && order_buffer[j].sec <sec_max){
      ready_to_scan=SCANBUF_SCANNING;
	  if((sec_max-order_buffer[j].sec)>1) printf("\nsec1: %d, sec_max: %d\n",order_buffer[j].sec,sec_max);
      shift_scale(j,order_buffer[j].ref_nr,order_buffer[j].ts0_ref);
      order_buffer[j].ts0_ref=0;
      order_buffer[j].flags=1;
    }
  }
  for(int j=0;j<MAXFEBNR+1;j++){  //scan seconds with max_sec
   if(order_buffer[j].flags==0){
     ready_to_scan=SCANBUF_SCANNING;
     shift_scale(j,order_buffer[j].ref_nr,order_buffer[j].ts0_ref);
     order_buffer[j].ts0_ref=0;
     order_buffer[j].flags=1;
   }
  }
  ready_to_fill=PROBUF_READY_TO_FILL;
}
//if a hole second is in the processing buffer the hits are copyed and scaled for scanning//////////////////////////////////////
void crt::pair_builder::shift_scale(int mac, int ref_nr, int ts0_ref){ //scale all events and store them in the scan buffer
  for(int i=0; i< number_ms[mac];i++){
    evbuf_scan[mac][i]=evbuf_scan[mac][ev_counter_scan[mac]-number_ms[mac]+i];
  }
  for(int i=0;i<ref_nr;i++){ //loop over all hits of a plane
   evbuf_scan[mac][i+number_ms[mac]]=evbuf_pro[mac][i];    //shift the hit from the pro to the scan buffer
   long scale0=((evbuf_pro[mac][i].ts0*1e9)/ts0_ref+0.5);  //calculate the scaling factor
   long scale1=((evbuf_pro[mac][i].ts1*1e9)/ts0_ref+0.5);  //calculate the scaling factor
   evbuf_scan[mac][i+number_ms[mac]].ts0_scaled=(int)scale0;    //store the scaled value
   evbuf_scan[mac][i+number_ms[mac]].ts1_scaled=(int)scale1;    //store the scaled value
   //evbuf_scan[mac][i+number_ms[mac]].ts0_scaled=evbuf_pro[mac][i].ts0;    //store the scaled value
   //evbuf_scan[mac][i+number_ms[mac]].ts1_scaled=evbuf_pro[mac][i].ts1;    //store the scaled value
   evbuf_scan[mac][i+number_ms[mac]].ts0_ref=ts0_ref; 
	  if(run_mode==RAW_MODE){
		  if(evbuf_scan[mac][i].mac5!=0xFFFF && (evbuf_scan[mac][i].flags!=5 && evbuf_scan[mac][i].flags!=7)){
			  int bufnr=0;
			  beam_ev[bufnr][0].mac5=evbuf_scan[mac][i].mac5;
			  beam_ev[bufnr][0].flags=evbuf_scan[mac][i].flags;
			  beam_ev[bufnr][0].ts0=evbuf_scan[mac][i].ts0;
			  beam_ev[bufnr][0].ts1=evbuf_scan[mac][i].ts1;
			  beam_ev[bufnr][0].lostcpu=evbuf_scan[mac][i].lostcpu;
			  beam_ev[bufnr][0].lostfpga=evbuf_scan[mac][i].lostfpga;
			  for(int amp=0; amp<32; amp++) beam_ev[bufnr][0].adc[amp]=evbuf_scan[mac][i].adc[amp];
			  beam_ev[bufnr][0].recover=evbuf_scan[mac][i].sec;
			  beam_ev[bufnr][0].nrtrigger=evbuf_scan[mac][i].sec;
			  beam_ev[bufnr][0].nrtrigger_11=evbuf_scan[mac][i].sec;

			  store_data(bufnr, 1);
			  ready_to_scan=SCANBUF_READY_TO_FILL;

		  }
		  //printf("test: mac: %d sec: %d\n",evbuf_scan[mac][i].mac5,evbuf_scan[mac][i].sec);
		  ready_to_scan=SCANBUF_READY_TO_FILL;
	  }
   //fprintf(text,"%d %d %d %d %d %d\n", evbuf_scan[mac][i].mac5, evbuf_scan[mac][i].flags,evbuf_scan[mac][i].ts0, evbuf_scan[mac][i].ts0_scaled, evbuf_scan[mac][i].adc[1],evbuf_scan[mac][i].sec );
 }
  // just for beam filtering...
	if(run_mode==FILTER_MODE){  //generate a special event with the additional information (e.g seconds)
  ts0_ref_event[0].mac5=evbuf_pro[mac][ref_nr].mac5;
  ts0_ref_event[0].flags=evbuf_pro[mac][ref_nr].flags;
  ts0_ref_event[0].lostcpu=evbuf_pro[mac][ref_nr].lostcpu;
  ts0_ref_event[0].lostfpga=evbuf_pro[mac][ref_nr].lostfpga;
  ts0_ref_event[0].ts0=evbuf_pro[mac][ref_nr].ts0;
  ts0_ref_event[0].ts1=evbuf_pro[mac][ref_nr].ts1;
  for(int amp=0; amp<32; amp++) ts0_ref_event[0].adc[amp]=evbuf_pro[mac][ref_nr].adc[amp];
  if(evbuf_pro[mac][ref_nr].flags!=7 && evbuf_pro[mac][ref_nr].flags!=5){
    printf("\n Error\n");
    for(int i=-2;i<5;i++){
      printf("%d %d %d %d %d %d\n", evbuf_pro[mac][ref_nr+i].mac5, evbuf_pro[mac][ref_nr+i].flags,evbuf_pro[mac][ref_nr+i].ts0, evbuf_pro[mac][ref_nr+i].ts1, evbuf_pro[mac][ref_nr+i].adc[0],evbuf_pro[mac][ref_nr+i].sec);
    }
    printf("Error\n");
  }
  ts0_ref_event[1].mac5=0xFFFF;
  ts0_ref_event[1].flags=1;
  ts0_ref_event[1].lostcpu=999;
  ts0_ref_event[1].lostfpga=0;
  ts0_ref_event[1].ts0=evbuf_pro[mac][ref_nr].sec;
  ts0_ref_event[1].ts1=evbuf_pro[mac][ref_nr].ms;
  for(int amp=0; amp<32; amp++) ts0_ref_event[1].adc[amp]=0;
	}
  //up to here...
  
 for(int i=ref_nr+1;i<ev_counter_mac[mac];i++){ //shift the rest to the beginnig of the pro buffer
  evbuf_pro[mac][i-ref_nr-1]=evbuf_pro[mac][i];
 }
 ev_counter_mac[mac]=ev_counter_mac[mac]-ref_nr-1;  //set the hit counter of the pro buffer
 ev_counter_scan[mac]=ref_nr+number_ms[mac]; //set the hit counter of the scan buffer
 if(run_mode==FILTER_MODE || run_mode==FILTER_PAIR_MODE || run_mode==TS1_CORR) { //filters first only beam events and then eventuall looking  for pairs
	 filter_buffer(mac);
	 scan_buffer_filter(mac); 
 }
	if(run_mode==PAIR_MODE) {//just looking for pairs
	 scan_buffer(mac); 
 }
}
void crt::pair_builder::scan_buffer_filter(int mac){  //scan over all events of one plane over one sec and search all possible coincidences
  //printf("scan feb nr %d\n", mac);
  last2_ms[mac]=0;
  
  uint32_t second = evbuf_scan[mac][ev_counter_scan[mac]-1].sec;
  if(ts1ref_second[second%100]!=second){
    ts1ref_counter[second%100]=0;
    ts1ref_second[second%100]=second;
  }
  if(verbose_!=0){
    if(previous_sec_scan[mac]==second) printf("\nequal seconds in %d: presecond: %d, this second: %d\n", mac,previous_sec_scan[mac], second);
    if(second-previous_sec_scan[mac]>1) printf("\nlost second in %d: presecond: %d, this second: %d\n", mac,previous_sec_scan[mac], second);
  }
  previous_sec_scan[mac]=second;
  int filter_counter=0;
  for(int i=last1_ms[mac];i<ev_counter_scan[mac];i++){
    //copy events in buffer for filtering
    evbuf_filter[mac][filter_counter]=evbuf_scan[mac][i];
    filter_counter++;
    
    //set last2_ms
    if(evbuf_scan[mac][i].ts0>=(1e9-2*MSOVERLAP) && i>number_ms[mac]){
      if(last2_ms[mac]==0){
        last2_ms[mac]=i;
        //printf("last2: %d \n", last2_ms[mac]);
      }
    }
    //set last1_ms
    if((evbuf_scan[mac][i].ts0>=(1e9-MSOVERLAP) && i>number_ms[mac])||i==(ev_counter_scan[mac]-1)){
      if(last2_ms[mac]==0) last2_ms[mac]=i;
      last1_ms[mac]=i;
      number_ms[mac]=ev_counter_scan[mac]-last2_ms[mac];
      last1_ms[mac]=last1_ms[mac]-last2_ms[mac];
      //printf("last1: %d \n", last1_ms[mac]);
      break;
    }
    if(evbuf_scan[mac][i].flags==10 || evbuf_scan[mac][i].flags==11){
       //printf("mac: %d, flag: %d, ts0: %d, ts1 %d sec: %d, tot events: %d\n", evbuf_scan[mac][i].mac5,evbuf_scan[mac][i].flags,evbuf_scan[mac][i].ts0,evbuf_scan[mac][i].ts1,evbuf_scan[mac][i].sec, ev_counter_scan[mac]);
      int new_ts1ref=1;
      for(int j=0;j<ts1ref_counter[second%100];j++){
        if(my_abs(evbuf_scan[mac][i].ts0,ts1ref_buffer[second%100][j])<20000) new_ts1ref++;
        //else printf("abs: %d ts1: %d second %d evsec: %d\n",abs(evbuf_scan[mac][i].ts1-ts1ref_buffer[second%100][j]),ts1ref_buffer[second%100][j], ts1ref_second[second%100],evbuf_scan[mac][i].sec);
      }
      if(new_ts1ref==1){
        //printf("new ts1: %d, number: %d\n", evbuf_scan[mac][i].ts0, ts1ref_counter[second%100]);
        ts1ref_buffer[second%100][ts1ref_counter[second%100]]=evbuf_scan[mac][i].ts0;
        ts1ref_counter[second%100]++;
      }
    }
    //printf("event.sec= %d   second= %d,  i=%d\n", evbuf_scan[mac][i].sec, second,i);
    //number_ms[mac]=ev_counter_scan[mac]-last2_ms[mac];   
  }
  ready_to_scan=SCANBUF_READY_TO_FILL;
  ev_counter_filter[mac]=filter_counter;
  //printf("ts1 counter: %d\n",ts1ref_counter[second%100]);
  //printf("\n");
  store_data_ts0_buffer(mac);
  //prozessed_hit_counter+=ev_counter_scan[mac]+1;
}
void crt::pair_builder::filter_buffer(int mac){
  //int MAXDIFFERENCE=(MAX_TIME_PASTBEAM-MAX_TIME_PREBEAM)/2;
  //int OFFSET=(MAX_TIME_PASTBEAM+MAX_TIME_PREBEAM)/2;
  int second = evbuf_filter[mac][ev_counter_filter[mac]-1].sec;
  int beam_ev_counter=0, ts1_ref_counter=0;
  int buf_to_send;
  //int ts1_ref_local[100];
  uint32_t ts1_ref_approved[100][2];
  int approved_ts1ref=0;
  for(int i=0; i<ev_counter_filter[mac];i++){
    for(int j=0;j<ts1ref_counter[second%100];j++){
      if(my_abs((evbuf_filter[mac][i].ts0_scaled),ts1ref_buffer[second%100][j])<200){
        //check if ts1_ref is also an event without flags in this second
        nrwm1[mac]++;
        ts1_ref_approved[approved_ts1ref][1]=0; //0= not recovered
        ts1_ref_approved[approved_ts1ref][0]=ts1ref_buffer[second%100][j];
        if(evbuf_filter[mac][i].flags!=11) {
          evbuf_filter[mac][i].flags = evbuf_filter[mac][i].flags | 0x100;
          ts1_ref_approved[approved_ts1ref][1]=1;}//1= recovered
        else nrwm2[mac]++;
        approved_ts1ref++;
      }
    }
  }
  if(run_mode==TS1_CORR){ //just assigne new ts1 for all events 6ms for beam trigger event.
    for(int i=0; i<ev_counter_filter[mac];i++){
      for(int j=0;j<approved_ts1ref;j++){
        if((ts1_ref_approved[j][0]-evbuf_filter[mac][i].ts0_scaled)<6e6 && (ts1_ref_approved[j][0]>evbuf_filter[mac][i].ts0_scaled)){     
            evbuf_filter[mac][i].ts1=4e9+ts1_ref_approved[j][0]-evbuf_filter[mac][i].ts0;
        }
      }
    }
    for(int i=0; i<ev_counter_filter[mac];i++){
      evbuf_filter_scan[mac][i]=evbuf_filter[mac][i];
      evbuf_filter_scan[mac][i].nrtrigger=nrwm2[mac];
      //evbuf_filter_scan[mac][i].nrtrigger_11=nrwm2[mac];
      //printf("1.i=%d mac %d, flag: %d, ts0 %d, adc[1] %d, sec: %d\n",i, evbuf_filter_scan[mac][i].mac5, evbuf_filter_scan[mac][i].flags,evbuf_filter_scan[mac][i].ts0, evbuf_filter_scan[mac][i].adc[1],evbuf_filter_scan[mac][i].sec);
    }
  }
  
  if(run_mode==FILTER_PAIR_MODE || run_mode==FILTER_MODE){
    for(int i=0; i<ev_counter_filter[mac];i++){
      for(int j=0;j<approved_ts1ref;j++){
        if(((ts1_ref_approved[j][0]>evbuf_filter[mac][i].ts0_scaled) && (ts1_ref_approved[j][0]-evbuf_filter[mac][i].ts0_scaled)<MAX_TIME_PREBEAM) ||
           ((ts1_ref_approved[j][0]<=evbuf_filter[mac][i].ts0_scaled) && (evbuf_filter[mac][i].ts0_scaled-ts1_ref_approved[j][0])<MAX_TIME_PASTBEAM)){

          evbuf_filter_scan[mac][beam_ev_counter]=evbuf_filter[mac][i];
          if(ts1_ref_approved[j][1]==1) evbuf_filter_scan[mac][beam_ev_counter].recover+=10;
          if((evbuf_filter[mac][i].ts0_scaled+1e9*(evbuf_filter[mac][i].sec-second))<ts1_ref_approved[j][0]){
            evbuf_filter_scan[mac][beam_ev_counter].ts1=4e9+ts1_ref_approved[j][0]-evbuf_filter[mac][i].ts0;}
          else evbuf_filter_scan[mac][beam_ev_counter].ts1=evbuf_filter[mac][i].ts1;
          evbuf_filter_scan[mac][beam_ev_counter].nrtrigger=nrwm1[mac];
          //evbuf_filter_scan[mac][beam_ev_counter].nrtrigger_11=ts0ref_counter[mac];
          beam_ev_counter++;
          break;
        }
      }
      if(evbuf_filter[mac][i].flags==11 || evbuf_filter[mac][i].flags==10) {
        //ts1_ref_local[ts1_ref_counter]=evbuf_filter[mac][i].ts0;
        ts1_ref_counter++;
      }
    }
  }
  if(run_mode==FILTER_PAIR_MODE) {
    ev_counter_filter_scan[mac]=beam_ev_counter;
	  scan_filter_buffer(mac); 
 }
  
  if(run_mode==TS1_CORR) { //scan all for pairs, with corrected ts1...
    ev_counter_filter_scan[mac]=ev_counter_filter[mac];
	  scan_filter_buffer(mac); 
  } 
  if(run_mode==FILTER_MODE){
    for(int i=0; i<beam_ev_counter; i++){
      beam_ev[send_bufnr][i].mac5=evbuf_filter_scan[mac][i].mac5;
      beam_ev[send_bufnr][i].flags=evbuf_filter_scan[mac][i].flags;
      beam_ev[send_bufnr][i].lostcpu=evbuf_filter_scan[mac][i].lostcpu;
      beam_ev[send_bufnr][i].lostfpga=evbuf_filter_scan[mac][i].lostfpga;
      beam_ev[send_bufnr][i].ts0=evbuf_filter_scan[mac][i].ts0;
      beam_ev[send_bufnr][i].ts1=evbuf_filter_scan[mac][i].ts1;
      for(int amp=0; amp<32; amp++) beam_ev[send_bufnr][i].adc[amp]=evbuf_filter_scan[mac][i].adc[amp];
      beam_ev[send_bufnr][i].recover=evbuf_filter_scan[mac][i].recover;
      beam_ev[send_bufnr][i].nrtrigger=evbuf_filter_scan[mac][i].nrtrigger;
      beam_ev[send_bufnr][i].nrtrigger_11=evbuf_filter_scan[mac][i].nrtrigger_11;
      //evbuf_filter_scan[mac][i]
      
    }
  beam_ev[send_bufnr][beam_ev_counter].mac5=0xFFFF;
  beam_ev[send_bufnr][beam_ev_counter].flags=beam_ev_counter;
  beam_ev[send_bufnr][beam_ev_counter].lostcpu=ts1_ref_counter;
  beam_ev[send_bufnr][beam_ev_counter].lostfpga=approved_ts1ref;
  beam_ev[send_bufnr][beam_ev_counter].ts0=evbuf_filter[mac][ev_counter_filter[mac]-1].sec;
  beam_ev[send_bufnr][beam_ev_counter].ts1=evbuf_filter[mac][ev_counter_filter[mac]-1].ms;

  for(int amp=0;amp<32;amp++) beam_ev[send_bufnr][beam_ev_counter].adc[amp]=0;
  beam_ev_counter++;
  buf_to_send=send_bufnr;
  //send_coinc(buf_to_send, j);
  store_data(buf_to_send, beam_ev_counter);
  store_data_ts0(mac);
  send_bufnr=(send_bufnr+1)%10;
  }
}
void crt::pair_builder::store_data_ts0_buffer(int mac){
  ts0_ref_event_buffer[mac][0].mac5=ts0_ref_event[0].mac5;
  ts0_ref_event_buffer[mac][0].flags=ts0_ref_event[0].flags;
  ts0_ref_event_buffer[mac][0].ts0=ts0_ref_event[0].ts0;
  ts0_ref_event_buffer[mac][0].ts1=ts0_ref_event[0].ts1;
  for(int j=0; j<32;j++) ts0_ref_event_buffer[mac][0].adc[j]=0;
  ts0_ref_event_buffer[mac][0].lostcpu=ts0_ref_event[0].lostcpu;
  ts0_ref_event_buffer[mac][0].lostfpga=ts0_ref_event[0].lostfpga;
  ts0_ref_event_buffer[mac][0].recover=0;
  ts0_ref_event_buffer[mac][0].nrtrigger=nrwm1[mac];
  ts0_ref_event_buffer[mac][0].nrtrigger_11=nrwm2[mac];
  //second event
  ts0_ref_event_buffer[mac][1].mac5=ts0_ref_event[1].mac5;
  ts0_ref_event_buffer[mac][1].flags=ts0_ref_event[1].flags;
  ts0_ref_event_buffer[mac][1].ts0=ts0_ref_event[1].ts0;
  ts0_ref_event_buffer[mac][1].ts1=ts0_ref_event[1].ts1;
  for(int j=0; j<32;j++) ts0_ref_event_buffer[mac][1].adc[j]=0;
  ts0_ref_event_buffer[mac][1].lostcpu=ts0_ref_event[1].lostcpu;
  ts0_ref_event_buffer[mac][1].lostfpga=ts0_ref_event[1].lostfpga;
  ts0_ref_event_buffer[mac][1].recover=0;
  ts0_ref_event_buffer[mac][1].nrtrigger=nrwm1[mac];
  ts0_ref_event_buffer[mac][1].nrtrigger_11=nrwm2[mac];
}
void crt::pair_builder::store_data_ts0(int mac){
  fwrite(ts0_ref_event_buffer[mac], sizeof(EVENT_t_send)*2, 1, file_store);  
}
void crt::pair_builder::store_data(int bufnr, int found_coinc){
  fwrite(beam_ev[bufnr], sizeof(EVENT_t_send)*found_coinc, 1, file_store);
  
}
//scanned all hits of one module and looks for time coincidences in the second of all other modules/////////////////////////////////////////
void crt::pair_builder::scan_buffer(int mac){  //scan over all events of one plane over one sec and search all possible coincidences
  //printf("scan feb nr %d\n", mac);
  if(evbuf_scan[mac][0].sec!=evbuf_scan[mac][ev_counter_scan[mac]-1].sec){
    for(int i=0; i<ev_counter_scan[mac]; i++) printf("mac %d, flag: %d, ts0: %d ts0_scaled %d, adc[1] %d, sec: %d, recover: %d\n", evbuf_scan[mac][i].mac5, evbuf_scan[mac][i].flags,evbuf_scan[mac][i].ts0, evbuf_scan[mac][i].ts1, evbuf_scan[mac][i].adc[1],evbuf_scan[mac][i].sec,evbuf_scan[mac][i].recover);
  }
  long time1, time2;
  long delta;
  //int coinc_counter=1;
  for(int i=0;i<ev_counter_scan[mac];i++){
    //coinc_counter=1;
    time1=evbuf_scan[mac][i].ts0_scaled;
    //printf("time1= %ld\n",time1);
    for(int j=0;j<MAXFEBNR;j++){
        if(j!=mac && ev_counter_scan[j]!=0){
        for(int k=0;k<ev_counter_scan[j];k++){
         //time2=(evbuf_scan[j][k].sec-evbuf_scan[mac][i].sec)*1e9+evbuf_scan[j][k].ts0_scaled;
         time2=evbuf_scan[j][k].ts0_scaled;
         //printf("time2= %ld\n",time2);
         delta=time2-time1;
          
         if((abs(delta)<MAX_TIME_DIFFERENCE)&&((evbuf_scan[j][k].flags==3 && evbuf_scan[mac][i].flags==3)||(evbuf_scan[j][k].flags==1 && evbuf_scan[mac][i].flags==1))){
           while(ready_to_send[send_bufnr]){ printf("wait bufnr: %d writting!\n", send_bufnr); }
           if(verbose_!=0){
             if(evbuf_scan[mac][i].sec<1496553028 && evbuf_scan[mac][i].sec>1496553021){
             printf("coincidence found: %ld delta\n",delta);
             printf("mac %d, flag: %d, ts0: %d ts0_scaled %d, adc[1] %d, sec: %d\n", evbuf_scan[mac][i].mac5, evbuf_scan[mac][i].flags,evbuf_scan[mac][i].ts0, evbuf_scan[mac][i].ts0_scaled, evbuf_scan[mac][i].adc[1],evbuf_scan[mac][i].sec);
             printf("mac %d, flag: %d, ts0 %d, ts0_scaled %d, adc[1] %d, sec: %d\n", evbuf_scan[j][k].mac5, evbuf_scan[j][k].flags,evbuf_scan[j][k].ts0, evbuf_scan[j][k].ts0_scaled, evbuf_scan[j][k].adc[1],evbuf_scan[j][k].sec);
             }
           }
           coincidence[send_bufnr][0].mac5=evbuf_scan[mac][i].mac5;
           coincidence[send_bufnr][0].flags=evbuf_scan[mac][i].flags;
           coincidence[send_bufnr][0].lostcpu=evbuf_scan[mac][i].lostcpu;
           coincidence[send_bufnr][0].lostfpga=evbuf_scan[mac][i].lostfpga;
           coincidence[send_bufnr][0].ts0=evbuf_scan[mac][i].ts0_scaled;
           coincidence[send_bufnr][0].ts1=evbuf_scan[mac][i].ts1;
           for(int amp=0;amp<32;amp++) coincidence[send_bufnr][0].adc[amp]=evbuf_scan[mac][i].adc[amp];
           coincidence[send_bufnr][0].recover=evbuf_scan[mac][i].recover;
           coincidence[send_bufnr][0].nrtrigger=evbuf_scan[mac][i].nrtrigger;
           coincidence[send_bufnr][0].nrtrigger_11=evbuf_scan[mac][i].nrtrigger_11;

           coincidence[send_bufnr][1].mac5=evbuf_scan[j][k].mac5;
           coincidence[send_bufnr][1].flags=evbuf_scan[j][k].flags;
           coincidence[send_bufnr][1].lostcpu=evbuf_scan[j][k].lostcpu;
           coincidence[send_bufnr][1].lostfpga=evbuf_scan[j][k].lostfpga;
           coincidence[send_bufnr][1].ts0=evbuf_scan[j][k].ts0_scaled;
           coincidence[send_bufnr][1].ts1=evbuf_scan[j][k].ts1;
           for(int amp=0;amp<32;amp++) coincidence[send_bufnr][1].adc[amp]=evbuf_scan[j][k].adc[amp];
           coincidence[send_bufnr][1].recover=evbuf_scan[j][k].recover;
           coincidence[send_bufnr][1].nrtrigger=evbuf_scan[j][k].nrtrigger;
           coincidence[send_bufnr][1].nrtrigger_11=evbuf_scan[j][k].nrtrigger_11; 
           
           coincidence[send_bufnr][2].mac5=0xFFFF;
           coincidence[send_bufnr][2].flags=2;
           coincidence[send_bufnr][2].ts0=evbuf_scan[mac][i].sec;
           coincidence[send_bufnr][2].ts1=evbuf_scan[j][k].sec;
           for(int amp=0;amp<32;amp++) coincidence[send_bufnr][2].adc[amp]=0;
           coincidence[send_bufnr][2].adc[0]=evbuf_scan[mac][i].ms;
           coincidence[send_bufnr][2].adc[1]=evbuf_scan[j][k].ms;
           coincidence[send_bufnr][2].recover=abs(i-j);
           
           int buf_to_send=send_bufnr;
           store_data_pairs(buf_to_send, 3);
         }
        } //end loop through all events of feb j
       }
      }//end loop over all febs
 }
  ready_to_scan=SCANBUF_READY_TO_FILL;
  //prozessed_hit_counter+=ev_counter_scan[mac]+1;
}

void crt::pair_builder::store_data_pairs(int bufnr, int found_coinc){
  fwrite(coincidence[bufnr], sizeof(EVENT_t_send)*found_coinc, 1, file_store);
  make2DHit(bufnr);
}

//scanned all hits of one module and looks for time coincidences in the second of all other modules/////////////////////////////////////////
void crt::pair_builder::scan_filter_buffer(int mac){  //scan over all events of one plane over one sec and search all possible coincidences
  //printf("scan feb nr %d\n", mac);
  long time1, time2;
  long delta;
  //int coinc_counter=1;
  for(int i=0;i<ev_counter_filter_scan[mac];i++){
    //printf("2.i=%d, mac %d, flag: %d, ts0 %d, adc[1] %d, sec: %d\n",i, evbuf_filter_scan[mac][i].mac5, evbuf_filter_scan[mac][i].flags,evbuf_filter_scan[mac][i].ts0, evbuf_filter_scan[mac][i].adc[1],evbuf_filter_scan[mac][i].sec);
    //coinc_counter=1;
    if(run_mode!=TS1_CORR || !(i>20 && evbuf_filter_scan[mac][i].ts0 > (1e9-MSOVERLAP))){  
    //if(1){
    time1=evbuf_filter_scan[mac][i].ts0_scaled;
    //printf("time1= %ld\n",time1);
    for(int j=0;j<MAXFEBNR;j++){
        if(j!=mac){
        for(int k=0;k<ev_counter_filter_scan[j];k++){
          if(run_mode!=TS1_CORR || !(k>20 && evbuf_filter_scan[j][k].ts0 > (1e9-MSOVERLAP))){  
           //time2=(evbuf_filter_scan[j][k].sec-evbuf_filter_scan[mac][i].sec)*1e9+evbuf_filter_scan[j][k].ts0_scaled;
           time2=evbuf_filter_scan[j][k].ts0_scaled;
           //printf("time2= %ld\n",time2);
           delta=time2-time1;

           if((abs(delta)<MAX_TIME_DIFFERENCE)&&((evbuf_filter_scan[j][k].flags==3 && evbuf_filter_scan[mac][i].flags==3)||(evbuf_filter_scan[j][k].flags==1 && evbuf_filter_scan[mac][i].flags==1))){
             while(ready_to_send[send_bufnr]){ printf("wait bufnr: %d writting!\n", send_bufnr); 
                                              //waiting_counter++;
                                             }
             //printf("coincidence found: %ld delta\n",delta);
             //printf("mac %d, flag: %d, ts0: %d ts0_scaled %d, adc[1] %d, sec: %d\n", evbuf_filter_scan[mac][i].mac5, evbuf_filter_scan[mac][i].flags,evbuf_filter_scan[mac][i].ts0, evbuf_filter_scan[mac][i].ts0_scaled, evbuf_filter_scan[mac][i].adc[1],evbuf_filter_scan[mac][i].sec);
             //printf("mac %d, flag: %d, ts0 %d, ts0_scaled %d, adc[1] %d, sec: %d\n", evbuf_filter_scan[j][k].mac5, evbuf_filter_scan[j][k].flags,evbuf_filter_scan[j][k].ts0, evbuf_filter_scan[j][k].ts0_scaled, evbuf_filter_scan[j][k].adc[1],evbuf_filter_scan[j][k].sec);
             coincidence[send_bufnr][0].mac5=evbuf_filter_scan[mac][i].mac5;
             coincidence[send_bufnr][0].flags=evbuf_filter_scan[mac][i].flags;
             coincidence[send_bufnr][0].lostcpu=evbuf_filter_scan[mac][i].lostcpu;
             coincidence[send_bufnr][0].lostfpga=evbuf_filter_scan[mac][i].lostfpga;
             coincidence[send_bufnr][0].ts0=evbuf_filter_scan[mac][i].ts0_scaled;
             coincidence[send_bufnr][0].ts1=evbuf_filter_scan[mac][i].ts1;
             for(int amp=0;amp<32;amp++) coincidence[send_bufnr][0].adc[amp]=evbuf_filter_scan[mac][i].adc[amp];
             coincidence[send_bufnr][0].recover=evbuf_filter_scan[mac][i].recover;
             coincidence[send_bufnr][0].nrtrigger=evbuf_filter_scan[mac][i].nrtrigger;
             coincidence[send_bufnr][0].nrtrigger_11=evbuf_filter_scan[mac][i].nrtrigger_11;

             coincidence[send_bufnr][1].mac5=evbuf_filter_scan[j][k].mac5;
             coincidence[send_bufnr][1].flags=evbuf_filter_scan[j][k].flags;
             coincidence[send_bufnr][1].lostcpu=evbuf_filter_scan[j][k].lostcpu;
             coincidence[send_bufnr][1].lostfpga=evbuf_filter_scan[j][k].lostfpga;
             coincidence[send_bufnr][1].ts0=evbuf_filter_scan[j][k].ts0_scaled;
             coincidence[send_bufnr][1].ts1=evbuf_filter_scan[j][k].ts1;
             for(int amp=0;amp<32;amp++) coincidence[send_bufnr][1].adc[amp]=evbuf_filter_scan[j][k].adc[amp];
             coincidence[send_bufnr][1].recover=evbuf_filter_scan[j][k].recover;
             coincidence[send_bufnr][1].nrtrigger=evbuf_filter_scan[j][k].nrtrigger;
             coincidence[send_bufnr][1].nrtrigger_11=evbuf_filter_scan[j][k].nrtrigger_11;

             coincidence[send_bufnr][2].mac5=0xFFFF;
             coincidence[send_bufnr][2].flags=2;
             coincidence[send_bufnr][2].ts0=evbuf_filter_scan[mac][i].sec;
             coincidence[send_bufnr][2].ts1=evbuf_filter_scan[j][k].sec;
             for(int amp=0;amp<32;amp++) coincidence[send_bufnr][2].adc[amp]=0;
             coincidence[send_bufnr][2].adc[0]=evbuf_filter_scan[mac][i].ms;
             coincidence[send_bufnr][2].adc[1]=evbuf_filter_scan[j][k].ms;
             coincidence[send_bufnr][2].recover=abs(i-j);
             int buf_to_send=send_bufnr;
             store_data_pairs(buf_to_send, 3);
           }
          }
        } //end loop through all events of feb j
       }
      }//end loop over all febs
 }
  }
  ready_to_scan=SCANBUF_READY_TO_FILL;
  //prozessed_hit_counter+=ev_counter_filter_scan[mac]+1;
}

unsigned int crt::pair_builder::my_abs(unsigned int a, unsigned int b){
  unsigned int c=0;
  if(a<b) c=b-a;
  else c=a-b;
  return c;
}


int crt::pair_builder::XYtracks(int order,int mac1, int mac2, int mac3, int mac4){
  switch(order){
    case 41: return XYtracks_Bottom_Top(mac1, mac2, mac3, mac4); break;
    case 14: return XYtracks_Bottom_Top(mac3, mac4, mac1, mac2); break;
    case 42: return XYtracks_Bottom_Pipe(mac1, mac2, mac3, mac4); break;
    case 24: return XYtracks_Bottom_Pipe(mac3, mac4, mac1, mac2); break;
    case 43: return XYtracks_Bottom_Feed(mac1, mac2, mac3, mac4); break;
    case 34: return XYtracks_Bottom_Feed(mac3, mac4, mac1, mac2); break;
    case 31: return XYtracks_Feed_Top(mac1, mac2, mac3, mac4); break;
    case 13: return XYtracks_Feed_Top(mac3, mac4, mac1, mac2); break;
    case 32: return XYtracks_Feed_Pipe(mac1, mac2, mac3, mac4); break;
    case 23: return XYtracks_Feed_Pipe(mac3, mac4, mac1, mac2); break;
    case 21: return XYtracks_Pipe_Top(mac1, mac2, mac3, mac4); break;
    case 12: return XYtracks_Pipe_Top(mac3, mac4, mac1, mac2); break;
  }
  return 0;
}
int crt::pair_builder::XY_pair(int pos, int mac1, int mac2){
  switch(pos){
    case 3: return XYtracks_Top(mac1,mac2); break;
    case 2: return XYtracks_Pipe(mac1,mac2); break;
    case 1: return XYtracks_Feed(mac1,mac2); break;
    case 0: return XYtracks_Bottom(mac1,mac2); break;
  }
 return 0;
  
}

int crt::pair_builder::XYtracks_Bottom(int mac3, int mac4){
  //2 event bottom, 1. event top
  //mac3=22;
  //mac4=11;
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  //int check2=0;
  if(mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==22 || mac3==23 ||mac3==24||mac3==12){
    if(mac4==11 || mac4==14 ||mac4==17||mac4==18||mac4==19){
      check1=1;}
  }
  if(mac3==11 || mac3==14 ||mac3==17||mac3==18||mac3==19){
    if(mac4==22 || mac4==23 ||mac4==24||mac4==12){
      check1=1;}
  } 
  if(check1==1) return 1;
  else {refused++; return 0;}
}
int crt::pair_builder::XYtracks_Top(int mac1, int mac2){
  //2 event bottom, 1. event top
  //mac3=22;
  //mac4=11;
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  //int check1=0;
  int check2=0;
  if(mac1==0||mac2==0) {refused++; return 0;}
  if(mac1==109 || mac1==105 ||mac1==106||mac1==107  ||mac1==125  ||mac1==126  ||mac1==128  ||mac1==111  ||mac1==112  ||mac1==108  ||mac1==124  ||mac1==123  ||mac1==122){
    if(mac2==113 || mac2==114 ||mac2==115||mac2==116  ||mac2==117  ||mac2==118  ||mac2==127  ||mac2==121  ||mac2==120  ||mac2==119  ||mac2==129){
      check2=1;}
  }
  if(mac1==113 || mac1==114 ||mac1==115||mac1==116  ||mac1==117  ||mac1==118  ||mac1==127  ||mac1==121  ||mac1==120  ||mac1==119  ||mac1==129){
    if(mac2==109 || mac2==105 ||mac2==106||mac2==107  ||mac2==125  ||mac2==126  ||mac2==128  ||mac2==111  ||mac2==112  ||mac2==108  ||mac2==124  ||mac2==123  ||mac2==122){
      check2=1;}
  }
  if(check2==1) return 1;
  else {refused++; return 0;}
}
int crt::pair_builder::XYtracks_Feed(int mac3, int mac4){
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  //int check2=0;
  
  if(mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==52||mac3==31||mac3==29||mac3==28||mac3==27||mac3==26||mac3==30){
    if(mac4==60||mac4==58||mac4==56||mac4==61||mac4==59||mac4==57){
      check1=1;}
  }
  if(mac3==60||mac3==58||mac3==56||mac3==61||mac3==59||mac3==57){
    if(mac4==52||mac4==31||mac4==29||mac4==28||mac4==27||mac4==26||mac4==30){
      check1=1;}
  }
  if(check1==1) return 1;
  else {refused++; return 0;} 
}
int crt::pair_builder::XYtracks_Pipe(int mac1, int mac2){
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  //int check1=0;
  int check2=0;
  
  if(mac1==0||mac2==0) {refused++; return 0;}
  if(mac1==37||mac1==33||mac1==34||mac1==35||mac1==36||mac1==32||mac1==38||mac1==39||mac1==40||mac1==41||mac1==42||mac1==43||mac1==44 ||mac1==45){
    if(mac2==53||mac2==54||mac2==55||mac2==15||mac2==16||mac2==20||mac2==21||mac2==46||mac2==47||mac2==48||mac2==49||mac2==50||mac2==51){
      check2=1;}
  }
  if(mac1==53||mac1==54||mac1==55||mac1==15||mac1==16||mac1==20||mac1==21||mac1==46||mac1==47||mac1==48||mac1==49||mac1==50||mac1==51){
    if(mac2==37||mac2==33||mac2==34||mac2==35||mac2==36||mac2==32||mac2==38||mac2==39||mac2==40||mac2==41||mac2==42||mac2==43||mac2==44 ||mac2==45){
      check2=1;}
  }
  if(check2==1) return 1;
  else {refused++; return 0;} 
}


int crt::pair_builder::XYtracks_Bottom_Top(int mac3, int mac4, int mac1, int mac2){
  //2 event bottom, 1. event top
  //mac3=22;
  //mac4=11;
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  int check2=0;
  if(mac1==0||mac2==0||mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==22 || mac3==23 ||mac3==24||mac3==12){
    if(mac4==11 || mac4==14 ||mac4==17||mac4==18||mac4==19){
      check1=1;}
  }
  if(mac3==11 || mac3==14 ||mac3==17||mac3==18||mac3==19){
    if(mac4==22 || mac4==23 ||mac4==24||mac4==12){
      check1=1;}
  } 
  if(mac1==109 || mac1==105 ||mac1==106||mac1==107  ||mac1==125  ||mac1==126  ||mac1==128  ||mac1==111  ||mac1==112  ||mac1==108  ||mac1==124  ||mac1==123  ||mac1==122){
    if(mac2==113 || mac2==114 ||mac2==115||mac2==116  ||mac2==117  ||mac2==118  ||mac2==127  ||mac2==121  ||mac2==120  ||mac2==119  ||mac2==129){
      check2=1;}
  }
  if(mac1==113 || mac1==114 ||mac1==115||mac1==116  ||mac1==117  ||mac1==118  ||mac1==127  ||mac1==121  ||mac1==120  ||mac1==119  ||mac1==129){
    if(mac2==109 || mac2==105 ||mac2==106||mac2==107  ||mac2==125  ||mac2==126  ||mac2==128  ||mac2==111  ||mac2==112  ||mac2==108  ||mac2==124  ||mac2==123  ||mac2==122){
      check2=1;}
  }
  if(check1==1 && check2==1) return 1;
  else {refused++; return 0;}
}

int crt::pair_builder::XYtracks_Bottom_Pipe(int mac3, int mac4, int mac1, int mac2){
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  int check2=0;
  if(mac1==0||mac2==0||mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==22 || mac3==23 ||mac3==24||mac3==12){
    if(mac4==11 || mac4==14 ||mac4==17||mac4==18||mac4==19){
      check1=1;}
  }
  if(mac3==11 || mac3==14 ||mac3==17||mac3==18||mac3==19){
    if(mac4==22 || mac4==23 ||mac4==24||mac4==12){
      check1=1;}
  } 
  if(mac1==37||mac1==33||mac1==34||mac1==35||mac1==36||mac1==32||mac1==38||mac1==39||mac1==40||mac1==41||mac1==42||mac1==43||mac1==44 ||mac1==45){
    if(mac2==53||mac2==54||mac2==55||mac2==15||mac2==16||mac2==20||mac2==21||mac2==46||mac2==47||mac2==48||mac2==49||mac2==50||mac2==51){
      check2=1;}
  }
  if(mac1==53||mac1==54||mac1==55||mac1==15||mac1==16||mac1==20||mac1==21||mac1==46||mac1==47||mac1==48||mac1==49||mac1==50||mac1==51){
    if(mac2==37||mac2==33||mac2==34||mac2==35||mac2==36||mac2==32||mac2==38||mac2==39||mac2==40||mac2==41||mac2==42||mac2==43||mac2==44 ||mac2==45){
      check2=1;}
  }
  if(check1==1 && check2==1) return 1;
  else {refused++; return 0;} 
}
int crt::pair_builder::XYtracks_Bottom_Feed(int mac3, int mac4, int mac1, int mac2){
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  int check2=0;
  if(mac1==0||mac2==0||mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==22 || mac3==23 ||mac3==24||mac3==12){
    if(mac4==11 || mac4==14 ||mac4==17||mac4==18||mac4==19){
      check1=1;}
  }
  if(mac3==11 || mac3==14 ||mac3==17||mac3==18||mac3==19){
    if(mac4==22 || mac4==23 ||mac4==24||mac4==12){
      check1=1;}
  } 
  if(mac1==52||mac1==31||mac1==29||mac1==28||mac1==27||mac1==26||mac1==30){
    if(mac2==60||mac2==58||mac2==56||mac2==61||mac2==59||mac2==57){
      check2=1;}
  }
  if(mac1==60||mac1==58||mac1==56||mac1==61||mac1==59||mac1==57){
    if(mac2==52||mac2==31||mac2==29||mac2==28||mac2==27||mac2==26||mac2==30){
      check2=1;}
  }
  if(check1==1 && check2==1) return 1;
  else {refused++; return 0;} 
}

int crt::pair_builder::XYtracks_Pipe_Top(int mac3, int mac4, int mac1, int mac2){
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  int check2=0;
  if(mac1==0||mac2==0||mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==37||mac3==33||mac3==34||mac3==35||mac3==36||mac3==32||mac3==38||mac3==39||mac3==40||mac3==41||mac3==42||mac3==43||mac3==44 ||mac3==45){
    if(mac4==53||mac4==54||mac4==55||mac4==15||mac4==16||mac4==20||mac4==21||mac4==46||mac4==47||mac4==48||mac4==49||mac4==50||mac4==51){
      check1=1;}
  }
  if(mac3==53||mac3==54||mac3==55||mac3==15||mac3==16||mac3==20||mac3==21||mac3==46||mac3==47||mac3==48||mac3==49||mac3==50||mac3==51){
    if(mac4==37||mac4==33||mac4==34||mac4==35||mac4==36||mac4==32||mac4==38||mac4==39||mac4==40||mac4==41||mac4==42||mac4==43||mac4==44 ||mac4==45){
      check1=1;}
  }
  if(mac1==109 || mac1==105 ||mac1==106||mac1==107  ||mac1==125  ||mac1==126  ||mac1==128  ||mac1==111  ||mac1==112  ||mac1==108  ||mac1==124  ||mac1==123  ||mac1==122){
    if(mac2==113 || mac2==114 ||mac2==115||mac2==116  ||mac2==117  ||mac2==118  ||mac2==127  ||mac2==121  ||mac2==120  ||mac2==119  ||mac2==129){
      check2=1;}
  }
  if(mac1==113 || mac1==114 ||mac1==115||mac1==116  ||mac1==117  ||mac1==118  ||mac1==127  ||mac1==121  ||mac1==120  ||mac1==119  ||mac1==129){
    if(mac2==109 || mac2==105 ||mac2==106||mac2==107  ||mac2==125  ||mac2==126  ||mac2==128  ||mac2==111  ||mac2==112  ||mac2==108  ||mac2==124  ||mac2==123  ||mac2==122){
      check2=1;}
  }
  if(check1==1 && check2==1) return 1;
  else {refused++; return 0;} 
}

int crt::pair_builder::XYtracks_Feed_Pipe(int mac3, int mac4, int mac1, int mac2){
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  int check2=0;
  
  if(mac1==0||mac2==0||mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==52||mac3==31||mac3==29||mac3==28||mac3==27||mac3==26||mac3==30){
    if(mac4==60||mac4==58||mac4==56||mac4==61||mac4==59||mac4==57){
      check1=1;}
  }
  if(mac3==60||mac3==58||mac3==56||mac3==61||mac3==59||mac3==57){
    if(mac4==52||mac4==31||mac4==29||mac4==28||mac4==27||mac4==26||mac4==30){
      check1=1;}
  }
  if(mac1==37||mac1==33||mac1==34||mac1==35||mac1==36||mac1==32||mac1==38||mac1==39||mac1==40||mac1==41||mac1==42||mac1==43||mac1==44 ||mac1==45){
    if(mac2==53||mac2==54||mac2==55||mac2==15||mac2==16||mac2==20||mac2==21||mac2==46||mac2==47||mac2==48||mac2==49||mac2==50||mac2==51){
      check2=1;}
  }
  if(mac1==53||mac1==54||mac1==55||mac1==15||mac1==16||mac1==20||mac1==21||mac1==46||mac1==47||mac1==48||mac1==49||mac1==50||mac1==51){
    if(mac2==37||mac2==33||mac2==34||mac2==35||mac2==36||mac2==32||mac2==38||mac2==39||mac2==40||mac2==41||mac2==42||mac2==43||mac2==44 ||mac2==45){
      check2=1;}
  }
  if(check1==1 && check2==1) return 1;
  else {refused++; return 0;} 
}

int crt::pair_builder::XYtracks_Feed_Top(int mac3, int mac4, int mac2, int mac1){
  //printf("mac1: %d mac2: %d mac3: %d mac4 %d\n", mac1,mac2,mac3,mac4);
  int check1=0;
  int check2=0;
  if(mac1==0||mac2==0||mac3==0||mac4==0) {refused++; return 0;}
  if(mac3==52||mac3==31||mac3==29||mac3==28||mac3==27||mac3==26||mac3==30){
    if(mac4==60||mac4==58||mac4==56||mac4==61||mac4==59||mac4==57){
      check1=1;}
  }
  if(mac3==60||mac3==58||mac3==56||mac3==61||mac3==59||mac3==57){
    if(mac4==52||mac4==31||mac4==29||mac4==28||mac4==27||mac4==26||mac4==30){
      check1=1;}
  }
  if(mac1==109 || mac1==105 ||mac1==106||mac1==107  ||mac1==125  ||mac1==126  ||mac1==128  ||mac1==111  ||mac1==112  ||mac1==108  ||mac1==124  ||mac1==123  ||mac1==122){
    if(mac2==113 || mac2==114 ||mac2==115||mac2==116  ||mac2==117  ||mac2==118  ||mac2==127  ||mac2==121  ||mac2==120  ||mac2==119  ||mac2==129){
      check2=1;}
  }
  if(mac1==113 || mac1==114 ||mac1==115||mac1==116  ||mac1==117  ||mac1==118  ||mac1==127  ||mac1==121  ||mac1==120  ||mac1==119  ||mac1==129){
    if(mac2==109 || mac2==105 ||mac2==106||mac2==107  ||mac2==125  ||mac2==126  ||mac2==128  ||mac2==111  ||mac2==112  ||mac2==108  ||mac2==124  ||mac2==123  ||mac2==122){
      check2=1;}
  }
  if(check1==1 && check2==1) return 1;
  else {refused++; return 0;} 
}


void crt::pair_builder::InitReconst(){
  printf("initializing the data!\n");
 double x, y, z, dt; //pos
 int id;//ind;
 int pl,ch,a,b,c;
  for(int m=0;m<200;m++) for(int c=0;c<32;c++) {Xs[m][c]=0;Ys[m][c]=0;Zs[m][c]=0;Dts[m]=0;Exists[m]=0;} 
  std::ifstream in;
  in.open("/uboone/app/users/dlorca/testberncode_june/larsoft_v06_36_00/srcs/uboonecode/uboone/CRT/CRTpositionsSiPM-V8.txt");
  if (in.is_open()){
    std::cout << "initializing position" << std::endl;
    while (!in.eof()) {
      in>>id>>x>>y>>z>>a>>b>>c;
      pl=id/100; ch=(id-100*int(pl));
     //ad-hoc correction of flipped modules for top
      //if(pl==108 || (pl>=111 && pl<=124) || pl==127|| pl==129) ch=31-ch; 
       //shift all to the detector coordinate system
      Xs[pl][ch]=x;
      Ys[pl][ch]=y;
      Zs[pl][ch]=z;
     // Exists[pl]=1;
    }      
    in.close();
  }
   else
  {
    std::cout << "Error opening position file";
  }
  
  std::ifstream in1;
  in1.open("/uboone/app/users/dlorca/testberncode_june/larsoft_v06_36_00/srcs/uboonecode/uboone/CRT/FEB_CableDelay-V8.txt");
  if (in1.is_open()){
    std::cout << "initializing delay" << std::endl;
    while (!in1.eof()) {
      in1>>pl>>dt;
      Dts[pl]=dt;
      Exists[pl]=1;
    }      
    in1.close();
    }
   else
   {
    std::cout << "Error opening cabledelay file";
   }
  
}
double getFineCoordExperim(int sL, int sR){
 double ret;
if(sL==0) sL=1;
ret=STRIPW/2.*atan(log(1.*sR/sL));

//if(ret>STRIPW/2) ret=STRIPW/2;
//if(ret<-STRIPW/2) ret=-STRIPW/2;
return ret;
}

double crt::pair_builder::TWCorrection(int max1ach, int max2ach) 
{
 return 4*log(max2ach/max1ach);
}

double crt::pair_builder::getHitT1(int mac1, int strip1,  int t1, int mac2, int strip2, int t2){ 
  double ret=0;
  double L=0;
  if(Xs[mac2][0]!=Xs[mac2][2]) L=Xs[mac2][strip2*2]-Xs[mac1][0]; //coord along the strip minus SiPM position
  if(Ys[mac2][0]!=Ys[mac2][2]) L=Ys[mac2][strip2*2]-Ys[mac1][0]; //coord along the strip minus SiPM position
  if(Zs[mac2][0]!=Zs[mac2][2]) L=Zs[mac2][strip2*2]-Zs[mac1][0]; //coord along the strip minus SiPM position
  if(L<0) L=-L;
  ret=t1-L*6.2/100.; //protagation along the fiber correction, 6.2 ns/m
  ret=ret+Dts[mac1]; //ref PPS cable 
  return ret;
}

double crt::pair_builder::getHitT2(int mac1, int strip1,  int t1, int mac2, int strip2, int t2){ 
  double ret=0; 
  double L=0;
  if(Xs[mac1][0]!=Xs[mac1][2]) L=Xs[mac1][strip1*2]-Xs[mac2][0]; //coord along the strip minus SiPM position
  if(Ys[mac1][0]!=Ys[mac1][2]) L=Ys[mac1][strip1*2]-Ys[mac2][0]; //coord along the strip minus SiPM position
  if(Zs[mac1][0]!=Zs[mac1][2]) L=Zs[mac1][strip1*2]-Zs[mac2][0]; //coord along the strip minus SiPM position
  if(L<0) L=-L;
  ret=t2-L*6.2/100.; //protagation along the fiber correction,6.2 ns/m
  ret=ret+Dts[mac2];//ref PPS cable  
  return ret;
}

double crt::pair_builder::getHitDT(int mac1, int strip1,  int t1, int max1_adc, int mac2, int strip2,int t2, int max2_adc){ 
 return getHitT2(mac1,strip1,t1,mac2,strip2,t2)-getHitT1(mac1,strip1,t1,mac2,strip2,t2)+4*log((double)max2_adc/8178)-4*log((double)max1_adc/8178);
}

double crt::pair_builder::getHitT(int mac1, int strip1,  int t1, int max1_adc, int mac2, int strip2,  int t2, int max2_adc){ 
 return (getHitT2(mac1,strip1,t1,mac2,strip2,t2)+getHitT1(mac1,strip1,t1,mac2,strip2,t2)+4*log((double)max1_adc/8178)+4*log((double)max2_adc/8178))/2.;
}

double crt::pair_builder::getHitX(int mac1, int strip1, int mac2, int strip2){ 
  double ret=0; 
  //printf("xs11: %f, Xs12: %f\n", Xs[mac1][0],Xs[mac1][2]);
  //printf("xs21: %f, Xs22: %f\n", Xs[mac1][0],Xs[mac1][2]);
  if(Xs[mac1][0]!=Xs[mac1][2]) ret=(Xs[mac1][strip1*2]+Xs[mac1][strip1*2+1])/2.;
  else if(Xs[mac2][0]!=Xs[mac2][2]) ret=(Xs[mac2][strip2*2]+Xs[mac2][strip2*2+1])/2.;
  else ret=(Xs[mac2][strip2*2]+Xs[mac1][strip1*2])/2;
  //printf("ret2: %f\n", ret);
  return ret;
}

double crt::pair_builder::getHitY(int mac1, int strip1, int mac2, int strip2){ 
  double ret=0; 
  if(Ys[mac1][0]!=Ys[mac1][2]) ret=(Ys[mac1][strip1*2]+Ys[mac1][strip1*2+1])/2.;
  else if(Ys[mac2][0]!=Ys[mac2][2]) ret=(Ys[mac2][strip2*2]+Ys[mac2][strip2*2+1])/2.;
  else ret=(Ys[mac2][strip2*2]+Ys[mac1][strip1*2])/2;
  return ret;
}

double crt::pair_builder::getHitZ(int mac1, int strip1, int mac2, int strip2){ 
  double ret=0; 
  if(Zs[mac1][0]!=Zs[mac1][2]) ret=(Zs[mac1][strip1*2]+Zs[mac1][strip1*2+1])/2.;
  else if(Zs[mac2][0]!=Zs[mac2][2]) ret=(Zs[mac2][strip2*2]+Zs[mac2][strip2*2+1])/2.;
  else ret=(Zs[mac2][strip2*2]+Zs[mac1][strip1*2])/2;
  return ret;
}

void crt::pair_builder::make2DHit(int bufnr){
  uint32_t first_second=0;
  crt::CRTHit crt2Dhit;
  if(coincidence[bufnr][2].ts0>1 && coincidence[bufnr][2].ts0<first_second) first_second=coincidence[bufnr][2].ts0;
  int XY = 0;
    crt2Dhit.plane=-1;
    if(XY_pair(0, coincidence[bufnr][0].mac5, coincidence[bufnr][1].mac5)){ XY=1; crt2Dhit.plane=0;}
    else if(XY_pair(1, coincidence[bufnr][0].mac5, coincidence[bufnr][1].mac5)){ XY=1; crt2Dhit.plane=1;}
    else if(XY_pair(2, coincidence[bufnr][0].mac5, coincidence[bufnr][1].mac5)){ XY=1; crt2Dhit.plane=2;}
    else if(XY_pair(3, coincidence[bufnr][0].mac5, coincidence[bufnr][1].mac5)){ XY=1; crt2Dhit.plane=3;}
  if(XY && coincidence[bufnr][2].ts0!=0){
    //if(XY){
    int max1_ach=0, max1_nch=0, max2_ach=0, max2_nch=0;  //max1_adc=0,  max2_adc=0

    crt2Dhit.feb_id.push_back(coincidence[bufnr][0].mac5);
    crt2Dhit.feb_id.push_back(coincidence[bufnr][1].mac5);

    std::vector<std::pair<int,double> > vec_pes_tevt;
    std::vector<std::pair<int,double> > vec_pes_st;

    for(size_t i_chan=0; i_chan<32; ++i_chan){
      /*int key_tevt = coincidence[bufnr][0].mac5*100+i_chan;
      std::pair<double,double> gain_tevt = crt::auxfunctions::getGain(key_tevt, SiPMgain);
      std::pair<double,double> pedestal_tevt = crt::auxfunctions::getGain(key_tevt, SiPMpedestal);
      double pes_tevt = ( (coincidence[bufnr][0].adc[i_chan]) - pedestal_tevt.first) / gain_tevt.first;*/
      double pes_tevt =  (double)coincidence[bufnr][0].adc[i_chan];
      std::pair<int,double> pair_tevt = std::make_pair(i_chan,pes_tevt);
      vec_pes_tevt.push_back(pair_tevt);

      /*int key_st = feb_st*100+i_chan;
      std::pair<double,double> gain_st = crt::auxfunctions::getGain(key_st, SiPMgain);
      std::pair<double,double> pedestal_st = crt::auxfunctions::getGain(key_st, SiPMpedestal);
      double pes_st = ( coincidence[bufnr][1].adc[i_chan] - pedestal_st.first) / gain_st.first;*/
      double pes_st =  (double)coincidence[bufnr][1].adc[i_chan];
      std::pair<int,double> pair_st = std::make_pair(i_chan,pes_st);
      vec_pes_st.push_back(pair_st);
    }


    std::map< uint8_t, std::vector<std::pair<int,double> > > pesmap_hit;
    pesmap_hit[coincidence[bufnr][0].mac5] = vec_pes_tevt;
    pesmap_hit[coincidence[bufnr][1].mac5] = vec_pes_st;
    crt2Dhit.pesmap = pesmap_hit;


    crt2Dhit.peshit=0;
    crt2Dhit.ts0_s=coincidence[bufnr][2].ts0;
    crt2Dhit.ts0_s_err=0;
    
    int ch[16];
    for(int j=0;j<16;j++){ //set max_ach + max__nch
      ch[j]=0;
      ch[j]=coincidence[bufnr][0].adc[j*2]+coincidence[bufnr][0].adc[j*2+1];
      if(ch[j]>max1_ach) {max1_ach=ch[j]; max1_nch=j;}
      ch[j]=0;
      ch[j]=coincidence[bufnr][1].adc[j*2]+coincidence[bufnr][1].adc[j*2+1];
      if(ch[j]>max2_ach) {max2_ach=ch[j]; max2_nch=j;}
    }
    int ts1_local1=0, ts1_local2=0;
    if(coincidence[bufnr][0].ts1>4e9){ ts1_local1=-(coincidence[bufnr][0].ts1-4e9);}
    else ts1_local1=coincidence[bufnr][0].ts1;
    if(coincidence[bufnr][1].ts1>4e9){ ts1_local2=-(coincidence[bufnr][1].ts1-4e9);}
    else ts1_local2=coincidence[bufnr][1].ts1;
    crt2Dhit.ts0_ns=getHitT(coincidence[bufnr][0].mac5, max1_nch,  coincidence[bufnr][0].ts0, max1_ach, coincidence[bufnr][1].mac5, max2_nch,  coincidence[bufnr][1].ts0, max2_ach);
    crt2Dhit.ts0_ns_err=abs((coincidence[bufnr][0].ts0+coincidence[bufnr][0].ts0)/2-crt2Dhit.ts0_ns);
    crt2Dhit.ts1_ns=getHitT(coincidence[bufnr][0].mac5, max1_nch,  ts1_local1, max1_ach, coincidence[bufnr][1].mac5, max2_nch,  ts1_local2, max2_ach);
    crt2Dhit.ts1_ns_err=abs((ts1_local1+ts1_local2)/2-crt2Dhit.ts1_ns);

    crt2Dhit.x_pos= getHitX(coincidence[bufnr][0].mac5, max1_nch, coincidence[bufnr][1].mac5, max2_nch);
    crt2Dhit.x_err=1;
    crt2Dhit.y_pos= getHitY(coincidence[bufnr][0].mac5, max1_nch, coincidence[bufnr][1].mac5, max2_nch);
    crt2Dhit.y_err=1;
    crt2Dhit.z_pos= getHitZ(coincidence[bufnr][0].mac5, max1_nch, coincidence[bufnr][1].mac5, max2_nch);
    crt2Dhit.z_err=1;

    //   printf("max1_ach: %d, max1_nch %d, max2_ach: %d, max2_nch: %d\n", max1_ach,max1_nch, max2_ach, max2_nch);
    //printf("ids: %d, %d, x: %f, y: %f, z: %f, plane: %d\n", crt2Dhit.feb_id[0], crt2Dhit.feb_id[1], crt2Dhit.x_pos, crt2Dhit.y_pos, crt2Dhit.z_pos, crt2Dhit.plane);
    //printf("ts0: %d, ts1: %d, sec: %d\n", crt2Dhit.ts0_s, crt2Dhit.ts1_ns, crt2Dhit.ts0_s);
  //something like this here:
    allmyCRTHits.push_back(crt2Dhit);
    //crt2Dhit.feb_id.clear();
  
}
}

