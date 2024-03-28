Kendall_cor<-function(dat1,dat2){
  c<-0
  d<-0
  t_x<-0
  t_y<-0
  for(i in 1:(length(dat1)-1)){
    for(j in (i+1):length(dat2)){
      if(((dat1[i]-dat1[j])*(dat2[i]-dat2[j]))>0){
        c=c+1
      }
      else if(((dat1[i]-dat1[j])*(dat2[i]-dat2[j]))<0){
        d=d+1
      }
      else if((dat1[i]-dat1[j])==0&(dat2[i]-dat2[j])!=0){
        t_x=t_x+1
      }
      else if((dat1[i]-dat1[j])!=0&(dat2[i]-dat2[j])==0){
        t_y=t_y+1
      }
    }
  }
  tau_b=(c-d)/sqrt((c+d+t_x)*(c+d+t_y))
  return(tau_b)
}
