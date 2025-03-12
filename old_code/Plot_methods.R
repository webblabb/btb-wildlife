#Categorical plotting function for bTB_wildlife model

#Usage of input parameters, defaults
#data_vec, totalPop_vec -- data values from the model. data vec is the category of interest while the other should usually be the N column from the output
                           #This is only used in scaling, and can use values other than N if desired. This is also optional if only doing regular plots. 
                           #otherwise this will throw an error

#herdSize, infType, dataType -- string arguments used in plot naming

#n_reps, -- number of replicates of the model in the data vectors

#n_years -- number of years run in the simulation

#minMax -- boolean, whether the min and max lines will be plotted. default False

#rep_Plots = 0 -- numerical argument, number of replicates of the model to be plotted. should be less than n_reps

#Scaled = F, Regular = T, Paired = F -- boolean, which types of plots will be created. Scaled gives the average proportion of the data_vect relative to totalPop.
                                       #Regular just plots the data vector mean over each time step. Paired generates both plots separately on the same plot window.
                                       #Defaults Regular = TRUE, Scaled, Paired = FALSE

#as.jpeg = F, to.console = T -- boolean, type of output from the function, writing the plot as a .jpeg or display in the R plot window. Both can be enabled simultaneously. 
                               #Defaults as.jpeg = FALSE, to.console = TRUE
#color -- designated R color name or value. Default blue -- may get more functionality in the future

writePlot <- function(data_vec = c(1), totalPop_vec = NULL, herdSize = "", infType = "", dataType = "", category = "", n_reps = 0, n_years = 0, minMax = F, rep_Plots = 0, Scaled = F, Regular = T, Paired = F, as.jpeg = F, to.console = T, color = 'blue'){
  
  #check for invalid input pairings
  #################################
  if(n_reps < rep_Plots)
  {
    print('Replicates to be plotted exceeds number of model replicates')
    reps = n_reps
  }
  
  if( (Scaled == T | Paired == T) & is.null(totalPop_vec) )
  {
    print('no value for total population vector (totalPop_vec) was given. Incompatible input with Scaled and Paired plot functionality')
    break
  }
  #################################
  
  
  
  #start method calculations for each plot type if called
  #######################################################
  
  #regular plots
  ##############
  if(Regular)
  {
    par(mfrow=c(1,1))
    data = data_vec
    n_tsteps = 12*n_years + 1
    t_vec = c(0:(n_tsteps-1)) #x variable for plots
    mean_vec = replicate( n_tsteps, 0 ) #storage vec for mean values at each time step
    lwr = 1 #lower bound replicate index
    upr = n_tsteps #upper bound index
    min_vec = data[ lwr:upr ] #use first replicate to ensure min comparisons are correct (i.e min would be all 0 if initialized to that)
    max_vec = data[ lwr:upr ] #shouldn't matter, but for the sake of consistency...
    
    for( i in 1:n_reps ) #loop over each replicate
    {
      lwr = 1 + (i-1)*n_tsteps #increment by replicate length each time except the first
      upr = i*n_tsteps #increment by number of replicates each time
      temp = data[ lwr:upr ]
      min_vec = pmin( min_vec, temp )
      max_vec = pmax( max_vec, temp )  
      
      for( t in 1:n_tsteps ) #loop over each time step
      {
        mean_vec[t] = mean_vec[t] + temp[t] #add the temp value to the vector of means -- more like a sum vector here
      }   
      
    }#end loop calculating mean, min and max
    
    mean_vec = mean_vec / n_reps #now take the mean at each time step with element-wise division
    y_min = min(min_vec) - abs(min(min_vec)*.05)
    y_max = max(max_vec) + max(max_vec)*.05
    
    
    ########################################################################################################
    # Write out plots in desired format ####################################################################
    ########################################################################################################
    ########################################################################################################
    #write to .jpeg loop
    ####################
    if(as.jpeg){
      fname = paste0()
      if(rep_Plots == 0){ #plot normally
        
        if(minMax){ #plot min max if enabled
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max))
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
        }else{ # normal plotting otherwise
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max)) 
        } 
        #end no rep plots conditional   
        
      }else{ #plot individual replicates
        
        if(minMax) #plots with everything
        {
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines(x = t_vec, y = temp, col='gray69', type='l') #add rep line
          }
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec, col=color, lwd=2 ) #ensure mean line is printed on top
          
        }else{ #plots with rep lines 
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample(x = c(1:n_reps), size = rep_Plots, replace = F)
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines(x = t_vec, y = temp, col='gray69', type='l') #add rep line
          }
          
          lines(x = t_vec, y = mean_vec, col=color, lwd=2) #ensure mean line is printed on top
        }
      }#end individual replicates loop
      dev.off()
    }
    ####################
    #end to jpeg loop
    
    
    #write to console loop
    ######################
    if(to.console){
      
      if(rep_Plots == 0){ #plot normally
        
        if(minMax){ #plot min max if enabled
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max))
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
        }else{ # normal plotting otherwise
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max)) 
        } 
      #end no rep plots conditional   
        
      }else{ #plot individual replicates
        
        if(minMax) #plots with everything
        {
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines( x = t_vec, y = temp, col='gray69', type='l' ) #add rep line
          }
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec, col=color, lwd=2 ) #ensure mean line is printed on top
          
        }else{ #plots with rep lines 
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines( x = t_vec, y = temp, col='gray69', type='l' ) #add rep line
          }
        
          lines( x = t_vec, y = mean_vec, col = color, lwd = 2 ) #ensure mean line is printed on top
        }
      }#end individual replicates loop
      
    }
    ######################  
    #end to.console loop  
    ########################################################################################################
  }#end regular plotting loop
  ##############
  
  
  #Scaled Plots
  #############
  if(Scaled) #note that variable names here are identical to normal plots, but the data variable is automatically scaled. this differs from plots with both
  {
    par(mfrow=c(1,1))
    data = data_vec / totalPop_vec
    data[is.nan(data)] <- 0
    n_tsteps = 12*n_years + 1
    t_vec = c(0:(n_tsteps-1)) #x variable for plots
    mean_vec = replicate( n_tsteps, 0 ) #storage vec for mean values at each time step
    lwr = 1 #lower bound replicate index
    upr = n_tsteps #upper bound index
    min_vec = data[ lwr:upr ] #use first replicate to ensure min comparisons are correct (i.e min would be all 0 if initialized to that)
    max_vec = data[ lwr:upr ] #shouldn't matter, but for the sake of consistency...
    
    for( i in 1:n_reps ) #loop over each replicate
    {
      lwr = 1 + (i-1)*n_tsteps #increment by replicate length each time except the first
      upr = i*n_tsteps #increment by number of replicates each time
      temp = data[ lwr:upr ]
      min_vec = pmin( min_vec, temp )
      max_vec = pmax( max_vec, temp )  
      
      for( t in 1:n_tsteps ) #loop over each time step
      {
        mean_vec[t] = mean_vec[t] + temp[t] #add the temp value to the vector of means -- more like a sum vector here
      }   
      
    }#end loop calculating mean, min and max
    
    mean_vec = mean_vec / n_reps #now take the mean at each time step with element-wise division
    y_min = -.05
    if(min(min_vec) < 0){y_min = min(min_vec) + min(min_vec)*.05}
    y_max = 1.05
    
    ########################################################################################################
    # Write out plots in desired format ####################################################################
    ########################################################################################################
    ########################################################################################################
    #write to .jpeg loop
    ####################
    if(as.jpeg){
      
      jpeg(paste0(fname,'.jpeg'))
      if(rep_Plots == 0){ #plot normally
        
        if(minMax){ #plot min max if enabled
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max))
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
        }else{ # normal plotting otherwise
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max)) 
        } 
        #end no rep plots conditional   
        
      }else{ #plot individual replicates
        
        if(minMax) #plots with everything
        {
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines(x = t_vec, y = temp, col='gray69', type='l') #add rep line
          }
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec, col=color, lwd=2 ) #ensure mean line is printed on top
          
        }else{ #plots with rep lines 
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample(x = c(1:n_reps), size = rep_Plots, replace = F)
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines(x = t_vec, y = temp, col='gray69', type='l') #add rep line
          }
          
          lines(x = t_vec, y = mean_vec, col=color, lwd=2) #ensure mean line is printed on top
        }
      }#end individual replicates loop
      dev.off()
    }
    ####################
    #end to jpeg loop
    
    
    #write to console loop
    ######################
    if(to.console){
      
      if(rep_Plots == 0){ #plot normally
        
        if(minMax){ #plot min max if enabled
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max))
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
        }else{ # normal plotting otherwise
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max)) 
        } 
        #end no rep plots conditional   
        
      }else{ #plot individual replicates
        
        if(minMax) #plots with everything
        {
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines( x = t_vec, y = temp, col='gray69', type='l' ) #add rep line
          }
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec, col=color, lwd=2 ) #ensure mean line is printed on top
          
        }else{ #plots with rep lines 
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines( x = t_vec, y = temp, col='gray69', type='l' ) #add rep line
          }
          
          lines( x = t_vec, y = mean_vec, col = color, lwd = 2 ) #ensure mean line is printed on top
        }
      }#end individual replicates loop
      
    }
    ######################  
    #end to.console loop  
    ########################################################################################################  
  }
  #############
  
  #Paired Plots
  #############
  if(Paired)
  {
    data = data_vec
    data_scaled = data_vec / totalPop_vec
    data_scaled[is.nan(data_scaled)] <- 0
    n_tsteps = 12*n_years + 1
    t_vec = c(0:(n_tsteps-1)) #x variable for plots
    mean_vec = replicate( n_tsteps, 0 ) #storage vec for mean values at each time step
    mean_vec_scaled = replicate( n_tsteps, 0 )
    lwr = 1 #lower bound replicate index
    upr = n_tsteps #upper bound index
    min_vec = data[ lwr:upr ] #use first replicate to ensure min comparisons are correct (i.e min would be all 0 if initialized to that)
    max_vec = data[ lwr:upr ] #shouldn't matter, but for the sake of consistency...
    min_vec_scaled = data_scaled[ lwr:upr ] #use first replicate to ensure min comparisons are correct (i.e min would be all 0 if initialized to that)
    max_vec_scaled = data_scaled[ lwr:upr ] #shouldn't matter, but for the sake of consistency...
    
    for( i in 1:n_reps ) #loop over each replicate
    {
      lwr = 1 + (i-1)*n_tsteps #increment by replicate length each time except the first
      upr = i*n_tsteps #increment by number of replicates each time
      
      temp = data[ lwr:upr ]
      min_vec = pmin( min_vec, temp )
      max_vec = pmax( max_vec, temp ) 
      
      temp_scaled = data_scaled[ lwr:upr ]
      min_vec_scaled = pmin( min_vec_scaled, temp_scaled )
      max_vec_scaled = pmax( max_vec_scaled, temp_scaled )  
      
      for( t in 1:n_tsteps ) #loop over each time step
      {
        mean_vec[t] = mean_vec[t] + temp[t] #add the temp value to the vector of means -- more like a sum vector here
        mean_vec_scaled[t] = mean_vec_scaled[t] + temp_scaled[t] #add the temp value to the vector of means -- more like a sum vector here
      }   
      
    }#end loop calculating mean, min and max
    
    mean_vec = mean_vec / n_reps #now take the mean at each time step with element-wise division
    y_min = min(min_vec) - abs(min(min_vec)*.05)
    y_max = max(max_vec) + .05*max(max_vec)
    
    mean_vec_scaled = mean_vec_scaled / n_reps #now take the mean at each time step with element-wise division
    y_min_scaled = -.05
    if(min(min_vec_scaled) < 0){y_min_scaled = min(min_vec_scaled) + min(min_vec_scaled)*.05}
    y_max_scaled = 1.05
    
    ########################################################################################################
    # Write out plots in desired format ####################################################################
    ########################################################################################################
    ########################################################################################################
    #write to .jpeg loop
    ####################
    if(as.jpeg)
    {
      
      jpeg(paste0(fname,'.jpeg'))
      par(mfrow=c(2,1))
      
      if(rep_Plots == 0){ #plot normally
        
        if(minMax){ #plot min max if enabled
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max))
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec_scaled, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min_scaled, y_max_scaled))
          lines( x = t_vec, y = min_vec_scaled, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec_scaled, col='black', type='l' ) #max line
          
        }else{ # normal plotting otherwise
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max)) 
          
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec_scaled, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min_scaled, y_max_scaled)) 
        } 
        #end no rep plots conditional   
        
      }else{ #plot individual replicates
        
        if(minMax) #plots with everything
        {
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines(x = t_vec, y = temp, col='gray69', type='l') #add rep line
          }
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec, col=color, lwd=2 ) #ensure mean line is printed on top
          
          
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec_scaled, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min_scaled, y_max_scaled))
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp_scaled=data_scaled[lwr:upr]
            lines(x = t_vec, y = temp_scaled, col='gray69', type='l') #add rep line
          }
          lines( x = t_vec, y = min_vec_scaled, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec_scaled, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec_scaled, col=color, lwd=2 ) #ensure mean line is printed on top
          
        }else{ #plots with rep lines 
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample(x = c(1:n_reps), size = rep_Plots, replace = F)
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines(x = t_vec, y = temp, col='gray69', type='l') #add rep line
          }
          lines(x = t_vec, y = mean_vec, col=color, lwd=2) #ensure mean line is printed on top
          
          
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec_scaled, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min_scaled, y_max_scaled))
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp_scaled=data_scaled[lwr:upr]
            lines(x = t_vec, y = temp_scaled, col='gray69', type='l') #add rep line
          }
          lines(x = t_vec, y = mean_vec_scaled, col=color, lwd=2) #ensure mean line is printed on top
          
        }
      }#end individual replicates loop
      dev.off()
    }
    ####################
    #end to jpeg loop
    
    
    #write to console loop
    ######################
    if(to.console){
      par(mfrow=c(2,1))
      if(rep_Plots == 0){ #plot normally
        
        if(minMax){ #plot min max if enabled
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max))
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec_scaled, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min_scaled, y_max_scaled))
          lines( x = t_vec, y = min_vec_scaled, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec_scaled, col='black', type='l' ) #max line
          
        }else{ # normal plotting otherwise
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min, y_max)) 
          
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot( x = t_vec, y = mean_vec_scaled, col=color, type = 'l', lwd = 2, main=plotName, xlab='time', ylab=category, ylim = c(y_min_scaled, y_max_scaled))
          
        } 
        #end no rep plots conditional   
        
      }else{ #plot individual replicates
        
        if(minMax) #plots with everything
        {
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines( x = t_vec, y = temp, col='gray69', type='l' ) #add rep line
          }
          lines( x = t_vec, y = min_vec, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec, col=color, lwd=2 ) #ensure mean line is printed on top
          
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec_scaled, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min_scaled, y_max_scaled))
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp_scaled=data_scaled[lwr:upr]
            lines( x = t_vec, y = temp_scaled, col='gray69', type='l' ) #add rep line
          }
          lines( x = t_vec, y = min_vec_scaled, col='black', type='l' ) #min line
          lines( x = t_vec, y = max_vec_scaled, col='black', type='l' ) #max line
          lines( x = t_vec, y = mean_vec_scaled, col=color, lwd=2 ) #ensure mean line is printed on top
          
        }else{ #plots with rep lines 
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min, y_max))
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          which_reps = sample( x = c(1:n_reps), size = rep_Plots, replace = F )
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp=data[lwr:upr]
            lines( x = t_vec, y = temp, col='gray69', type='l' ) #add rep line
          }
          
          lines( x = t_vec, y = mean_vec, col = color, lwd = 2 ) #ensure mean line is printed on top
          
          fname = paste0()
          plotName = paste0('Herd size: ',herdSize, ' - scaled ', dataType, ' with ', infType)
          plot(x = t_vec, y = mean_vec_scaled, col = color, type = 'l', lwd = 2, main = plotName, xlab = 'time', ylab = category, ylim = c(y_min_scaled, y_max_scaled))
          lwr = 1 #lower bound replicate index
          upr = n_tsteps #upper bound index
          for(i in 1:rep_Plots)
          {
            lwr = 1 + ((which_reps[i])-1)*n_tsteps #increment by replicate length each time except the first
            upr = (which_reps[i])*n_tsteps #increment by number of replicates each time
            temp_scaled=data_scaled[lwr:upr]
            lines( x = t_vec, y = temp_scaled, col='gray69', type='l' ) #add rep line
          }
          
          lines( x = t_vec, y = mean_vec_scaled, col = color, lwd = 2 ) #ensure mean line is printed on top
        }
      }#end individual replicates loop
      
    }
    ######################  
    #end to.console loop  
    ########################################################################################################
  }
  #############
  #######################################################
}
  
  #ask me why I put in so many toggles
tot_births=run$`S Births` + run$`SS_S Births`
norm_deaths = run$`S Deaths` + run$`E1 Deaths` + run$`I Deaths`
tot_deaths=run$`S Deaths` + run$`E1 Deaths` + run$`I Deaths` + run$`SS_S Deaths` + run$`SS_E1 Deaths` + run$`SS_I Deaths`
tot_inf = run$SS_I + run$I 

writePlot(data_vec = tot_inf, totalPop_vec = run$N, herdSize = "", infType = "", dataType = "", category = "SS_frac", n_reps = 50, n_years = 10, minMax = T, rep_Plots = 5, Scaled = T, Regular = T, Paired = T, as.jpeg = F, to.console = T, color = 'blue')
data_vec = run$`S Births`; totalPop_vec = tot_births; herdSize = " 500 "; infType = " inf "; dataType = " dat "; n_reps = 50; n_years = 10; minMax = T; rep_Plots = 5; Scaled = T; Regular = T; Paired = T; as.jpeg = F; to.console = T; color = 'blue'; category = "SS_frac"



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#plot function to plot several mean simulation results at once i.e. results over multiple categories.
multiPlots <- function(data = data.frame(), data_scaled = data.frame(), pairwiseComp = '', pairwiseVec = c(0), category = '', infType = '', n_reps = 0,  n_years = 0, scaled = F, regular = T, paired = F, col_vec = NA ){
  library(RColorBrewer)
  
  if( ( ncol(data) != ncol(data_scaled) ) && (paired || scaled) ){print('column numbers do not match'); break}
  if( ( ncol(data) > 12 ) && is.na(col_vec) ){print('too many data columns. colors must be manually assigned'); break}
  
  n_cols = ncol(data)
  colors = col_vec
  if(is.na(col_vec)){ colors=brewer.pal(n_cols,'Paired') }
  
  
  ########################################################################################
  ########################################################################################
  
  if(regular){
    
    par(mfrow=c(1,1))
    mean_dat = data.frame() #storage for full data means
    n_tsteps = 12*n_years + 1
    t_vec = c(0:(n_tsteps-1)) #x variable for plots
    mean_vec = replicate( n_tsteps, 0 ) #storage vec for single simulation mean values at each time step
    lwr = 1 #lower bound replicate index
    upr = n_tsteps #upper bound index
    min_vec = data[lwr:upr,1] #min value at each time step over all simulations - for plot range
    max_vec= data[lwr:upr,1] #max value at each time step over all simulations - for plot range
    
    for(i in 1:n_cols)
    {
      
      for( j in 1:n_reps ) #loop over each replicate within the ith simulation
      {
        lwr = 1 + (j-1)*n_tsteps #increment by replicate length each time except the first
        upr = j*n_tsteps #increment by number of replicates each time
        temp = data[ lwr:upr , i] #single replicate of simulation i data
        min_vec = pmin( min_vec, temp )
        max_vec = pmax( max_vec, temp ) 
      
        # for( t in 1:n_tsteps ) #loop over each time step
        # {
        #   mean_vec[t] = mean_vec[t] + temp[t] #add the temp value to the vector of means -- more like a sum vector here
        # }  
        #since vector addition is elementwise, there was no need for the loop
        mean_vec = mean_vec + temp #add the temp value to the vector of means -- more like a sum vector here
        
        
      }#end loop calculating mean
    
    mean_vec = mean_vec / n_reps #now take the mean at each time step with element-wise division
    y_min = min(min_vec) - abs(min(min_vec)*.05)
    y_max = max(max_vec) + abs(max(max_vec)*.05)
    
    if(i == 1){mean_dat <- as.data.frame(mean_vec)} # workaround to inability to merge vector with empty data frame - allows for unspecified size
    mean_dat[,i] <- mean_vec #assign mean values for simulation i to the storage dataframe
    mean_vec = replicate( n_tsteps, 0 ) #reset the mean vector for next simulation set
  
    }#end loop over columns
    
    ##start plot generation##
    plotName = paste0('Mean ' ,category , ' over several ', pairwiseComp)
    
    y_max = max(mean_dat) + .05*abs(max(mean_dat))
    plot( x = t_vec, y = mean_dat[,1], xaxt = "n",  col = colors[1], type = 'l', lwd = 2, main=plotName, xlab='time', ylab = category, xlim = c( 0, 1.2*max(t_vec) ), ylim = c(y_min, y_max) ) 
    axis(1, at = seq(0, max(t_vec), by = 12) )
    for(k in 2:n_cols){
      lines( x = t_vec, y = mean_dat[,k], col = colors[k], type='l' )
    }
    legend('topright', title = pairwiseComp, box.lty=0, inset=c(.03,.01), legend = pairwiseVec, fill=colors)
    
  }#end regular plot loop
  
  ########################################################################################
  ########################################################################################
  
  if(scaled){
    
    par(mfrow=c(1,1))
    mean_dat_scaled = data.frame() #storage for full data means
    n_tsteps = 12*n_years + 1
    t_vec = c(0:(n_tsteps-1)) #x variable for plots
    mean_vec_scaled = replicate( n_tsteps, 0 ) #storage vec for single simulation mean values at each time step
    lwr = 1 #lower bound replicate index
    upr = n_tsteps #upper bound index
    min_vec = data_scaled[lwr:upr,1] #min value at each time step over all simulations - for plot range
    max_vec= data_scaled[lwr:upr,1] #max value at each time step over all simulations - for plot range
    
    for(i in 1:n_cols)
    {
      
      for( j in 1:n_reps ) #loop over each replicate within the ith simulation
      {
        lwr = 1 + (j-1)*n_tsteps #increment by replicate length each time except the first
        upr = j*n_tsteps #increment by number of replicates each time
        temp = data_scaled[ lwr:upr , i] #single replicate of simulation i data
        min_vec_scaled = pmin( min_vec_scaled, temp )
        max_vec_scaled = pmax( max_vec_scaled, temp ) 
        
        # for( t in 1:n_tsteps ) #loop over each time step
        # {
        #   mean_vec_scaled[t] = mean_vec_scaled[t] + temp[t] #add the temp value to the vector of means -- more like a sum vector here
        # }   
        #since vector addition is elementwise, there was no need for the loop
        mean_vec = mean_vec + temp #add the temp value to the vector of means -- more like a sum vector here
        
        
      }#end loop calculating mean
      
      mean_vec = mean_vec / n_reps #now take the mean at each time step with element-wise division
      y_min_scaled = min(min_vec_scaled) - abs(min(min_vec_scaled)*.05)
      y_max_scaled = max(max_vec_scaled) + abs(max(max_vec_scaled)*.05)
      
      if(i == 1){mean_dat_scaled <- as.data.frame(mean_vec_scaled)} # workaround to inability to merge vector with empty data frame - allows for unspecified size
      mean_dat_scaled[,i] <- mean_vec_scaled #assign mean values for simulation i to the storage dataframe
      mean_vec_scaled = replicate( n_tsteps, 0 ) #reset the mean vector for next simulation set
      
    }#end loop over columns
    
    ##start plot generation##
    plotName = paste0('Mean scaled' ,category , ' over several ', pairwiseComp)
    y_max = max(max(mean_dat_scaled) + .05*abs(max(mean_dat_scaled)), 1.05)
    plot( x = t_vec, y = mean_dat_scaled[,1], xaxt = "n",  col = colors[1], type = 'l', lwd = 2, main=plotName, xlab='time', ylab = category, xlim = c( 0, 1.2*max(t_vec) ), ylim = c( min(0,y_min), y_max ) ) 
    axis(1, at = seq(1, max(t_vec), by = 12) )
    for(k in 2:n_cols){
      lines( x = t_vec, y = mean_dat_scaled[,k], col = colors[k], type='l' )
    }
    legend('topright', title = pairwiseComp, box.lty=0, inset=c(.03,.01), legend = pairwiseVec, fill=colors)
    
    
  }#end scaled plot loop
  
  ########################################################################################
  ########################################################################################
  
  if(paired){
    
    par(mfrow=c(2,1))
    
    mean_dat = data.frame() #storage for full data means
    mean_dat_scaled = data.frame() #storage for full data means
    
    n_tsteps = 12*n_years + 1
    t_vec = c(0:(n_tsteps-1)) #x variable for plots
    
    mean_vec = replicate( n_tsteps, 0 ) #storage vec for single simulation mean values at each time step
    mean_vec_scaled = replicate( n_tsteps, 0 ) #storage vec for single simulation mean values at each time step
    
    lwr = 1 #lower bound replicate index
    upr = n_tsteps #upper bound index
    
    min_vec = pmin( min_vec, temp )
    max_vec = pmax( max_vec, temp )
    min_vec_scaled = pmin( min_vec_scaled, temp )
    max_vec_scaled = pmax( max_vec_scaled, temp )
    
    for(i in 1:n_cols)
    {
      
      for( j in 1:n_reps ) #loop over each replicate within the ith simulation
      {
        lwr = 1 + (j-1)*n_tsteps #increment by replicate length each time except the first
        upr = j*n_tsteps #increment by number of replicates each time
        
        temp = data[ lwr:upr , i] #single replicate of simulation i data
        min_vec = pmin( min_vec, temp )
        max_vec = pmax( max_vec, temp ) 
        
        temp = data_scaled[ lwr:upr , i] #single replicate of simulation i data
        min_vec_scaled = pmin( min_vec_scaled, temp )
        max_vec_scaled = pmax( max_vec_scaled, temp ) 
        
        
        # for( t in 1:n_tsteps ) #loop over each time step
        # {
        #   mean_vec_scaled[t] = mean_vec_scaled[t] + temp[t] #add the temp value to the vector of means -- more like a sum vector here
        #   
        #   mean_vec[t] = mean_vec[t] + temp[t] #add the temp value to the vector of means -- more like a sum vector here
        # } 
        #since vector addition is elementwise, there was no need for the loop
        mean_vec_scaled = mean_vec_scaled + temp #add the temp value to the vector of means -- more like a sum vector here
        mean_vec = mean_vec + temp #add the temp value to the vector of means -- more like a sum vector here
        
        
      }#end loop calculating mean
      
      mean_vec = mean_vec / n_reps #now take the mean at each time step with element-wise division
      y_min = min(min_vec) - abs(min(min_vec)*.05)
      y_max = max(max_vec) + abs(max(max_vec)*.05)
      
      mean_vec_scaled = mean_vec_scaled / n_reps #now take the mean at each time step with element-wise division
      y_min_scaled = min(min_vec_scaled) - abs(min(min_vec_scaled)*.05)
      y_max_scaled = max(max_vec_scaled) + abs(max(max_vec_scaled)*.05)
      
      if(i == 1){mean_dat_scaled <- as.data.frame(mean_vec_scaled); mean_dat <- as.data.frame(mean_vec)} # workaround to inability to merge vector with empty data frame - allows for unspecified size
      
      mean_dat[,i] <- mean_vec #assign mean values for simulation i to the storage dataframe
      mean_vec = replicate( n_tsteps, 0 ) #reset the mean vector for next simulation set
      
      mean_dat_scaled[,i] <- mean_vec_scaled #assign mean values for simulation i to the storage dataframe
      mean_vec_scaled = replicate( n_tsteps, 0 ) #reset the mean vector for next simulation set
      
      
    }#end loop over columns
    
    ##start plot generation##
    par(mfrow=c(1,1))
    plotName = paste0('Mean ' ,category , ' over several ', pairwiseComp)
    y_max = max(mean_dat) + .05*abs(max(mean_dat))
    plot( x = t_vec, y = mean_dat[,1], xaxt = "n",  col = colors[1], type = 'l', lwd = 2, main=plotName, xlab='time', ylab = category, xlim = c( 0, 1.2*max(t_vec) ), ylim = c(y_min, y_max) ) 
    axis(1, at = seq(0, max(t_vec), by = 12) )
    for(k in 2:n_cols){
      lines( x = t_vec, y = mean_dat[,k], col = colors[k], type='l' )
    }
    legend('topright', title = pairwiseComp, box.lty=0, inset=c(.03,.01), legend = pairwiseVec, fill=colors)
    
    
    plotName = paste0('Mean scaled' ,category , ' over several ', pairwiseComp, 's')
    y_max_scaled = max(max(mean_dat_scaled) + .05*abs(max(mean_dat_scaled)), 1.05)
    plot( x = t_vec, y = mean_dat_scaled[,1], xaxt = "n",  col = colors[1], type = 'l', lwd = 2, main=plotName, xlab='time', ylab = category, xlim = c( 0, 1.2*max(t_vec) ), ylim = c( min(0,y_min_scaled), y_max_scaled) ) 
    axis(1, at = seq(0, max(t_vec), by = 12) )
    for(k in 2:n_cols){
      lines( x = t_vec, y = mean_dat_scaled[,k], col = colors[k], type='l' )
    }
    legend('topright', title = pairwiseComp, box.lty=0, inset=c(.03,.01), legend = pairwiseVec, fill=colors)
    
  }#end paired plot loop
    
    
  
  
  
}#end function

#data for this function needs to be input as a list of n data frames where n is the number of herd sizes.
#each data frame should be composed of m columns where m is the number of different initial infected individuals. these values must also be passed 
#to the function in the n_inf parameter. The dataframe columns should be the count of all infected individuals over all replicates of the simulation.
#in the future, there will be a function that neatly arranges the data properly to input to this function. 
fadeoutPlots <- function(data = list(), regular = T, fadeTime = F, paired = F, infType = '', n_reps = 0,  n_years = 0, n_inf = c(0), herdSizes = c(0), col_vec = NA){
  library(RColorBrewer)
  
  if( ( length(data) > 12 ) && is.na(colors) ){print('too many data columns. colors must be manually assigned'); break}
  
  
  
  
  ########################################################################################
  ########################################################################################
  fadeoutProp = data.frame() #n x m data frame for mean fadeout values
  fadeoutTime = data.frame() #n x m data frame for mean fadeout values
  
  for(l in 1:length(data)){#loop over the different herd sizes
    #initialize simulation characteristics and all storage variables
    #some will need to be set outside of the herd size
    InfData = data[[l]] #data frame from list - results for all number initial infected at a fixed herd size
    n_cols = ncol(InfData)
    fadeout_vec = replicate(length(n_inf),0)
    fadeTimes = replicate(length(n_inf),0)
    n_tsteps = 12*n_years + 1
    lwr = 1 #lower bound replicate index
    upr = n_tsteps #upper bound index
    reInf=0 #not currently in use
    for(i in 1:n_cols) #loop over simulations with changing number of initial infected
    {
    
      for( j in 1:n_reps ) #loop over each replicate within the ith simulation
      {
        lwr = 1 + (j-1)*n_tsteps #increment by replicate length each time except the first
        upr = j*n_tsteps #increment by number of replicates each time
        temp = InfData[ lwr:upr , i] #single replicate of simulation i data
        infTime = min(which(temp != 0))
      
        #apologies for the convoluted conditionals - this was necessary to ensure accurate fadeTime calculations
        if( sum(temp==0) >= 1){ #see if the number of infected reaches 0
          if( sum( (which(temp==0) > infTime) == T) > 0 ){ #only a fadeout if there was actually an infection before reaching 0 infections
            if( sum(temp[min(which(temp==0)[which(temp==0) > infTime]):length(temp)]) > 0){reInf = reInf+1} #test for reinfection after first fadeout event #not currently in use
            fadeout_vec[i] = fadeout_vec[i] + 1 #sum values over each replicate
            fadeTimes[i] = fadeTimes[i] + min(which(temp==0)[which(temp==0) > infTime]) #sum values over each replicate
          }
        }
      }#end loop over reps
      fadeout_vec[i] = fadeout_vec[i] / n_reps #compute means from iterative sums
      fadeTimes[i] = fadeTimes[i] / n_reps #compute means from iterative sums
    
    }#end loop over column
    if(l==1){fadeoutTime <- data.frame(fadeTimes); fadeoutProp <- data.frame(fadeout_vec)}
    fadeoutTime[,l] <- fadeTimes
    fadeoutProp[,l] <- fadeout_vec
    #fadeout_vec = replicate(length(n_inf),0) #reset for next herd size
    #fadeTimes = replicate(length(n_inf),0) #reset for next herd size
  }#end loop over list elements i.e. herd sizes
  
  #plots
  #############################################
  #############################################
  colors = col_vec
  if(is.na(col_vec)){ colors=brewer.pal(n_cols,'Paired') }
  
  if(regular){
    par(mfrow=c(1,1))
    plot(x=n_inf, y=fadeoutProp[,1], type='p', col=colors[1], xaxt = "n",xlim=c(1,1.2*max(n_inf)), ylim=c(0,min(1, 1.05*max(fadeoutProp))), pch = 17, main='Fadeout Probability by Number of Initial Infections', xlab='No. initial infected', ylab='prop. fadeout')
    axis( 1, at = n_inf )
    legend('topright', title = 'pop. size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
    for(l in 2:length(data)){
      lines(x=n_inf, y=fadeoutProp[,l], type='p', col=colors[l], pch = 17)
    }
  }#end regular plots
  
  
  if(fadeTime){
    par(mfrow=c(1,1))
    plot(x=n_inf, y=fadeoutTime[,1], type='p', col=colors[1], xaxt = "n", xlim=c(1,1.2*max(n_inf)), ylim=c(0, min(12*n_years, 1.05*max(fadeoutTime) )), pch = 17, main='Mean Fadeout time by Number of Initial Infections', xlab='No. initial infected', ylab='months')
    axis( 1, at = n_inf )
    legend('topright', title = 'pop. size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
    for(l in 2:length(data)){
      lines(x=n_inf, y=fadeoutTime[,l], type='p', col=colors[l], pch = 17)
    }
  }#end fadeTime plots
  
  
  if(paired){
    par(mfrow=c(2,1))
    
    plot(x=n_inf, y=fadeoutProp[,1], type='p', col=colors[1], xaxt = "n", xlim=c(1,1.2*max(n_inf)), ylim=c(0,min(1, 1.05*max(fadeoutProp))), pch = 17, main='Fadeout Probability by Number of Initial Infections', xlab='No. initial infected', ylab='prop. fadeout')
    axis( 1, at = n_inf )
    legend('topright', title = 'pop. size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
    for(l in 2:length(data)){
      lines(x=n_inf, y=fadeoutProp[,l], type='p', col=colors[l], pch = 17)
    }
    
    plot(x=n_inf, y=fadeoutTime[,1], type='p', col=colors[1], xaxt = "n", xlim=c(1,1.2*max(n_inf)), ylim=c(0, min(12*n_years, 1.05*max(fadeoutTime) )), pch = 17, main='Mean Fadeout time by Number of Initial Infections', xlab='No. initial infected', ylab='months')
    axis( 1, at = n_inf )
    legend('topright', title = 'pop. size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
    for(l in 2:length(data)){
      lines(x=n_inf, y=fadeoutTime[,l], type='p', col=colors[l], pch = 17)
    }
  }#end paired plots

  
}#end function
  



#data for this function is formatted as a list of data frames. each list element indicates a different herd size. 
#The dataframe contained in each element should contain n columns where n is the number of hunt mortality values used.
#these columns should contain scaled model output over the whole simulation of: 
#prop. infected hunted[t] / prop. infected total[t]
#an additional column for the month may also be included, since hunting only occurs for 3 months out of the year.
#as such, much of the data will be 0 or inf if there are no infections, and this is corrected for within the function
#paired data is similarly formatted but contains data on number of total infections
hunt_prop_plots <- function(data = list(), pairedData = list(NA), regular = T, paired = F, infType = '', n_reps = 0,  n_years = 0, par_vals = c(0), herdSizes = c(0), col_vec = NA){
  library(RColorBrewer)
  
  if( ( length(data) > 12 ) && is.na(colors) ){print('too many data columns. colors must be manually assigned'); break}
  if(is.na(pairedData[1]) && paired){print('Error: no paired data provided'); break}
  n_cols = length(data)
  colors = col_vec
  if(is.na(col_vec)){ colors=brewer.pal(n_cols,'Paired') }
  
  detection_Data=data.frame()
  ##################################################
  ##################################################
  for(l in 1:length(data))
  {
    
    huntData = data[[l]]
    n_cols = ncol(huntData)
    n_tsteps = 12*n_years + 1
    detectionRatio = replicate(n_cols ,0) #initialize or reset detection ratio for each different herd size
    
    for(i in 1:n_cols) #loop over simulations with changing number of initial infected
    {
      
      
      for( j in 1:n_reps ) #loop over each replicate within the ith simulation
      {
        lwr = 1 + (j-1)*n_tsteps #increment by replicate length each time except the first
        upr = j*n_tsteps #increment by number of replicates each time
        temp = huntData[ lwr:upr , i] #single replicate of simulation i data
        if(sum(temp %in% c(Inf,NA) == T) > 0 ){print('Warning: NA or Inf value detected')}
        temp = temp[which( !(temp %in% c(0,Inf,NA,NaN)) )] #remove undesired values that may occur i.e. not hunting season (0) or no infections (NaN). inf and NA entries are removed just in case
        detectionRatio[i] = detectionRatio[i] + sum(temp)/length(temp)
        
        
      }#end loop over reps
      detectionRatio = detectionRatio / n_reps #compute mean from iterative sums
      if(i == 1){detection_Data[,l] = detectionRatio}
      detection_Data[,l] = detectionRatio
      
    }#end loop over columns
    
  }#end loop. over herd sizes
  
  ##################################################
  ##################################################
  
  ##Plots##
  if(regular){
    par(mfrow=c(1,1))
    plot(x=par_vals, y=detection_Data[,1], type='p', col=colors[1], xaxt = "n",xlim=c(0,1.2*max(par_vals)), ylim=c(0,max(detection_Data)*1.05), pch = 17, main='detection ratio by hunting rate', xlab='hunting rate', ylab='ratio of observed:true prevalence')
    axis( 1, at = par_vals )
    lines(x=par_vals, y=replicate(length(par_vals),1), type = 'l', col = 'green') #y=1 line - i.e. ideal case where observed prop and true prop are equal.
    legend('topright', title = 'Herd Size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
    for(l in 2:length(data)){
      lines(x=par_vals, y=detection_Data[,l], type='p', col=colors[l], pch = 17)
    }
  }
  
}  



#function to reformat necessary run data for plots
#parameter type should be a length 3 vector with entries: 
#('multi' 'fadeout' or 'hunt') , ('seeded.csv' 'spillover.csv' 'seeded_SSinf.csv' or 'spillover_SSinf.csv') , (data column name from output)
#scaling is an option specifically for multi plots that will return a scaled version of the data
#by is a string giving the column name to scale by.
#path is a string argument with the path to the folder containing the run data to be reformatted.
#herdSizes is a vector of herd sizes to group by. this is a necessary argument for fadeout and hunting plots

#warning messages are expected, as dataframes are added to an empty list.
plotsData_reformat <- function(version = 2, type = c('','seeded.csv',''), runPath = '', scaled = F, scaleBy = NULL, herdSizes = c(0)){
  
  if( scaled == T && is.null(scaleBy) ){Print('Error: scaled enabled with no default "scaleBy" value'); break}
  
  
  #format data for type multi
  if(type[1] %in% c('multi', 'Multi')){
    outputData = data.frame() #create storage data for output
    file_list = intersect(list.files(path = runPath, pattern = type[2]), list.files(path = runPath, pattern = paste0('bTBwl_v.',version) ) )
    
    
    for(i in 1:length(file_list)){
      tempDat = read.csv(file = paste0(runPath, file_list[i])) # read ith run data
      tempVec = tempDat[,which(colnames(tempDat) %in% type[3])]
      if(scaled){
        scaleVec= tempDat[,which(colnames(tempDat) %in% scaleBy)]
        tempVec = tempVec / scaleVec
      }
      if(i==1){outputData = data.frame(tempVec)} #assign proper length to data on first iteration
      
      outputData[,i] = tempVec #save ith run data to dataframe
      
    }#end loop over files
    
  }#end multi data loop
  
  
  
  #format data for type fadeout w/ number infection
  if(type[1] %in% c('fadeout', 'Fadeout')){

    outputData = vector('list', length = length(herdSizes)) #create list for data frame storage
    
    for(k in 1:length(herdSizes)){
      tempData = data.frame() #create storage data for output
      file_list = intersect( list.files( path = runPath, pattern = paste0('_N', herdSizes[k], '_') ), 
                             intersect(list.files( path = runPath , pattern = paste0('#inf_', type[2]) ) , 
                                       list.files(path = runPath, pattern = paste0('bTBwl_v.',version) ) 
                                       ) 
                             )#file list with desired run type and particular herd size
     
      
      for(i in 1:length(file_list)){
        runDat = read.csv(file = paste0(runPath, file_list[i]) ) #read in run data. in this case this loop will iterate ofver different initial infected numbers for a fixed herd size
        
        tempVec = runDat$E1 + runDat$I + runDat$SS_E1 + runDat$SS_I #add up total number of infecteds over all reps
        if(i == 1){tempData = data.frame(tempVec)} #recreate data frame with temp vector length assigned
        tempData[,i] = tempVec #populate columns of tempData with infection counts for fixed herd size with varied number of infected
      }#end loop over individual runs
      tempData = as.data.frame(tempData)
      outputData[[k]] = as.data.frame(tempData)
    }#end loop over herd sizes
    
  }#end fadeout data loop
  
  #format data for type fadeout w/ percent infection
  if(type[1] %in% c('fadeout%', 'Fadeout%')){
    
    outputData = vector('list', length = length(herdSizes)) #create list for data frame storage
    
    for(k in 1:length(herdSizes)){
      tempData = data.frame() #create storage data for output
      file_list = intersect( list.files( path = runPath, pattern = paste0('_N', herdSizes[k], '_') ), 
                             intersect(list.files( path = runPath , pattern = paste0('%inf_', type[2]) ) ,
                                       list.files(path = runPath, pattern = paste0('bTBwl_v.',version) ) 
                                       ) 
                             ) #file list with desired run type and particular herd size
      
      
      for(i in 1:length(file_list)){
        runDat = read.csv(file = paste0(runPath, file_list[i]) ) #read in run data. in this case this loop will iterate ofver different initial infected numbers for a fixed herd size
        
        tempVec = runDat$E1 + runDat$I + runDat$SS_E1 + runDat$SS_I #add up total number of infecteds over all reps
        if(i == 1){tempData = data.frame(tempVec)} #recreate data frame with temp vector length assigned
        tempData[,i] = tempVec #populate columns of tempData with infection counts for fixed herd size with varied number of infected
      }#end loop over individual runs
      tempData = as.data.frame(tempData)
      outputData[[k]] = as.data.frame(tempData)
    }#end loop over herd sizes
    
  }#end fadeout data loop
  
  
  
  #format data for type hunt
  if(type[1] %in% c('hunt', 'Hunt')){
    
    outputData = vector('list', length = length(herdSizes)) #create list for data frame storage
    
    for(k in 1:length(herdSizes)){
      tempData = data.frame() #create storage data for output
      file_list =  intersect( list.files( path = runPath , pattern = type[2] ) , 
                              list.files(path = runPath, pattern = paste0('bTBwl_v.',version) ) 
                              ) #file list with desired run type and particular herd size
      for(i in 1:length(file_list)){
        runDat = read.csv(file = paste0(runPath, file_list[i])) #read in run data. in this case this loop will iterate ofver different initial infected numbers for a fixed herd size
        
        tot_hunted = runDat$E1.hunted + runDat$SS_E1.hunted + runDat$`I.hunted` + runDat$`SS_I.hunted` + runDat$`S.hunted` + runDat$`SS_S.hunted`
        infHunted = runDat$`E1.hunted` + runDat$`SS_E1.hunted` + runDat$`I.hunted` + runDat$`SS_I.hunted`
        totInf = runDat$E1 + runDat$I + runDat$SS_E1 + runDat$SS_I
        
        obsRatio =  infHunted / tot_hunted
        trueRatio = totInf / runDat$N
        
        
        tempVec = obsRatio / trueRatio #ratio of observed and true proportions
        
        if(i == 1){tempData = data.frame(tempVec)} #recreate data frame with temp vector length assigned
        tempData[,i] = tempVec #populate columns of tempData with infection counts for fixed herd size with varied number of infected
        
      }#end loop over individual runs
      tempData = as.data.frame(tempData)
      outputData[[k]] = as.data.frame(tempData)
    }#end loop over herd sizes
    
  }#end hunt data loop
  
  
  return(outputData) 
  
}#end function 


################################################################################
################################################################################
################################################################################

#formatting for hunting rate plots - for different herd sizes or hunting rates

##reformat run data into lists. - list by herd sizes
setwd('~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_huntRates/')
dat_by_huntRates = vector('list', length = length(hunt_rates))

for(i in 1:length(hunt_rates)){
  n_tsteps = years*3
  file_list = intersect(list.files(pattern = paste0('_', N[i], '_huntRate')) , 
                        intersect(list.files(pattern = paste0( 'seeded_2%inf.csv')),
                                  list.files( pattern = paste0('bTBwl_v.',vsn) )
                                  )
                        )
  #file_list = intersect(list.files(pattern = paste0('_', N[i], '_huntRate')) , list.files(pattern = paste0( 'seeded.csv')) )
  
  for(j in 1:length(file_list)){
    data = read.csv(file = file_list[j])
    data = data[which(data$quarter %in% 4),] #take only data where hunting can occur
    data = data[,c('N', 'E1', 'SS_E1', 'I', 'SS_I', 'E1.hunted', 'SS_E1.hunted', 'I.hunted', 'SS_I.hunted', 'Total_hunt', 'Total_inf')]
    
    if(j==1){out_data = as.data.frame(data$N[1:n_tsteps])} #assign something to the dataframe in first loop iteration
    
    #define true and observed prevalence
    Ntot = data$N
    tru_prop = (data$I + data$E1 + data$SS_I + data$SS_E1) / Ntot
    obs_prop = (data$SS_E1.hunted + data$SS_E1.hunted + data$E1.hunted + data$I.hunted) / data$Total_hunt
    
    #remove NaNs 0/0
    tru_prop[is.nan(tru_prop)] = 0
    obs_prop[is.nan(obs_prop)] = 0
    
    #now average over all simulations
    n_tsteps = years*3
    grandMean = c(0,0,0)
    avg_N = replicate(n_tsteps, 0)
    avg_tru_pop = replicate(n_tsteps, 0)
    avg_obs_prev = replicate(n_tsteps, 0)
    scaled_steps = 3*years
    for( n in 1:n_reps ) #loop over each replicate within the ith simulation
    {
      lwr = 1 + (n-1)*scaled_steps #increment by replicate length each time except the first
      upr = n*scaled_steps #increment by number of replicates each time
      temp1 = Ntot[ lwr:upr ] #single replicate of simulation i data
      temp2 = tru_prop[ lwr:upr ]
      temp3 = obs_prop[ lwr:upr ]
      
      avg_N =  avg_N + temp1
      avg_tru_pop = avg_tru_pop + temp2
      avg_obs_prev = avg_obs_prev + temp3
    }
    
    avg_N =  avg_N / n_reps
    avg_tru_pop = avg_tru_pop / n_reps
    avg_obs_prev = avg_obs_prev / n_reps
    
    #assign data and col names
    out_data[,3*j-2] = avg_N
    out_data[,3*j-1] = avg_tru_pop
    out_data[,3*j] = avg_obs_prev
    names(out_data)[3*j-2] = paste0('N_',j)
    names(out_data)[3*j-1] = paste0('true_prop_',j)
    names(out_data)[3*j] = paste0('obs_prop_',j)
  }
  #assign data frame to the list
  dat_by_herdSize[[i]] = out_data
  
}


##reformat run data into lists. - list by hunting rate
setwd('~/Documents/Brandon/bTB_wildlife_code/bTB_wildlife_runs/Var_huntRates/')
dat_by_huntRates = vector('list', length = length(hunt_rates))

for(i in 1:length(hunt_rates)){
  n_tsteps = years*3
  #file_list = intersect(list.files(pattern = paste0('_huntRate', hunt_rates[i],'_')) , list.files(pattern = paste0( 'seeded.csv')) )
  file_list = intersect(list.files(pattern = paste0('_huntRate', hunt_rates[i],'_')) , 
                        intersect(list.files(pattern = paste0( 'seeded_2%inf.csv')) ,
                                  list.files( pattern = paste0('bTBwl_v.',vsn) )
                                  )
                        )
  
  for(j in 1:length(file_list)){
    data = read.csv(file = file_list[j])
    data = data[which(data$quarter %in% 4),] #take only data where hunting can occur
    data = data[,c('N', 'E1', 'SS_E1', 'I', 'SS_I', 'E1.hunted', 'SS_E1.hunted', 'I.hunted', 'SS_I.hunted', 'Total_hunt', 'Total_inf')]
    
    if(j==1){out_data = as.data.frame(data$N[1:n_tsteps])} #assign something to the dataframe in first loop iteration
    
    #define true and observed prevalence
    Ntot = data$N
    tru_prop = (data$I + data$E1 + data$SS_I + data$SS_E1) / Ntot
    obs_prop = (data$SS_E1.hunted + data$SS_E1.hunted + data$E1.hunted + data$I.hunted) / data$Total_hunt
    
    #remove NaNs 0/0
    tru_prop[is.nan(tru_prop)] = 0
    obs_prop[is.nan(obs_prop)] = 0
    
    #now average over all simulations
    
    grandMean = c(0,0,0)
    avg_N = replicate(n_tsteps, 0)
    avg_tru_pop = replicate(n_tsteps, 0)
    avg_obs_prev = replicate(n_tsteps, 0)
    scaled_steps = 3*years
    for( n in 1:n_reps ) #loop over each replicate within the ith simulation
    {
      lwr = 1 + (n-1)*(scaled_steps) #increment by replicate length each time except the first
      upr = n*scaled_steps #increment by number of replicates each time
      temp1 = Ntot[ lwr:upr ] #single replicate of simulation i data
      temp2 = tru_prop[ lwr:upr ]
      temp3 = obs_prop[ lwr:upr ]
      
      avg_N =  avg_N + temp1
      avg_tru_pop = avg_tru_pop + temp2
      avg_obs_prev = avg_obs_prev + temp3
    }
    
    avg_N =  avg_N / n_reps
    avg_tru_pop = avg_tru_pop / n_reps
    avg_obs_prev = avg_obs_prev / n_reps
    
    #assign data and col names
    out_data[,3*j-2] = avg_N
    out_data[,3*j-1] = avg_tru_pop
    out_data[,3*j] = avg_obs_prev
    
    names(out_data)[3*j-2] = paste0('N_',j)
    names(out_data)[3*j-1] = paste0('true_prop_',j)
    names(out_data)[3*j] = paste0('obs_prop_',j)
  
  }
  #assign data frame to the list
  dat_by_huntRates[[i]] = out_data
  
}




#plots using above data
#by herd size...


simLength = 3*years 
t_vec = c(1:nrow(dat_by_herdSize[[1]]))
par(mfrow=c(1,1))
sizesToPlot = c(1,3,5,7)
sizesToPlot = c(1:7)
dat = dat_by_huntRates[[7]]
colors = brewer.pal(length(sizesToPlot),'Paired')
plot(x= t_vec, y=dat[,2], type='l', col=colors[1], xaxt = "n",xlim=c(0,1.2*max(t_vec)),  ylim = c(0,1), main='', xlab='', ylab='')
lines(x= t_vec, y=dat[,3], type='p', col=colors[1],  pch = 17)
axis( 1, at = t_vec )
#legend('topright', title = 'Herd Size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
for(i in 2:length(sizesToPlot)){
  
  lines(x= t_vec, y=dat[,3*i-1], type='l', col=colors[i]  ) 
  lines(x= t_vec, y=dat[,3*i], type='p', col=colors[i],  pch = 17)

}

par(mfrow=c(1,1))
ratesToPlot = c(1,3,5,7)
ratesToPlot = c(1:7)
dat = dat_by_herdSize[[5]]
colors = brewer.pal(length(ratesToPlot),'Paired')
plot(x= t_vec, y=dat[,2], type='l', col=colors[1], xaxt = "n", lwd = 2  ,xlim=c(0,1.2*max(t_vec)),  ylim = c(0,1), main='', xlab='', ylab='')
lines(x= t_vec, y=dat[,3], type='p', col=colors[1],  pch = 17)
axis( 1, at = t_vec )
#legend('topright', title = 'Herd Size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
for(i in 2:length(ratesToPlot)){
  
  lines(x= t_vec, y=dat[,3*i-1], type='l', col=colors[i], lwd = 2  ) 
  lines(x= t_vec, y=dat[,3*i], type='p', col=colors[i],  pch = 17)
  
}


#using ratio of observed to true proportion
for(w in 1:length(dat_by_herdSize)){
  dat = dat_by_herdSize[[w]]
  ratesToPlot = c(1:7)

  colors=ghibli::ghibli_palettes$PonyoMedium
  plot(x= t_vec, y= dat[,3] / dat[,2], type='l', col=colors[1], xaxt = "n", lwd = 2  ,xlim=c(0,1.2*max(t_vec)),  ylim = c(0,1.3), main=paste0('Prevalence estimates with population size ', herdSizes[w]), xlab='time', ylab='observed:true prevalence')
  
  axis( 1, at = t_vec )
  abline(a=1, b=0, col='black', lwd=2,lty='dashed', xlim = c(0,max(t_vec)), )
  legend('topright', title = 'Hunt Rate', box.lty = 0, inset = c(.03,.01), legend = hunt_rates, fill = colors)
  for(i in 2:length(ratesToPlot)){
  
    lines(x= t_vec, y=dat[,3*i]/dat[,3*i-1], type='l', col=colors[i], lwd = 2  ) 
 
  }
  
  plot(x= t_vec, y=replicate(length(t_vec), median(dat[,3])/median(dat[,2])), type='l', col=colors[1], xaxt = "n", lwd = 2  ,xlim=c(0,1.2*max(t_vec)),  ylim = c(0,1.3), main=paste0('Median prevalence estimate with population size ', herdSizes[w]), xlab='time', ylab='observed:true prevalence')
  
  axis( 1, at = t_vec )
  abline(a=1, b=0, col='black', lwd=2,lty='dashed', xlim = c(0,max(t_vec)), )
  legend('topright', title = 'Hunt Rate', box.lty = 0, inset = c(.03,.01), legend = hunt_rates, fill = colors)
  for(i in 2:length(ratesToPlot)){
    
    lines(x= t_vec, y=replicate(length(t_vec), median(dat[,3*i])/median(dat[,3*i-1])), type='l', col=colors[i], lwd = 2  ) 
    
  }
  
}


#using ratio of observed to true proportion
for(w in 1:length(dat_by_huntRates)){
  dat = dat_by_huntRates[[w]]
  ratesToPlot = c(1:7)
  
  colors=ghibli::ghibli_palettes$PonyoMedium
  yvec=dat[,3] / dat[,2]
  yvec[which(is.nan(yvec))] = 0
  plot(x= t_vec, y= dat[,3] / dat[,2], type='l', col=colors[1], xaxt = "n", lwd = 2  ,xlim=c(0,1.2*max(t_vec)),  ylim = c(0,1.3), main=paste0("Prevalence estimates with a hunting rate of ", hunt_rates[w]), xlab='time', ylab='observed:true prevalence')
  
  axis( 1, at = t_vec )
  abline(a=1, b=0, col='black', lwd=2, lty='dashed', xlim = c(0,10))
  legend('topright', title = 'Herd Size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
  for(i in 2:length(ratesToPlot)){
    
    lines(x= t_vec, y=dat[,3*i]/dat[,3*i-1], type='l', col=colors[i], lwd = 2  ) 
   
  }
  
  
  plot(x= t_vec, y= replicate(length(t_vec), median(dat[,3])/median(dat[,2])), type='l', col=colors[1], xaxt = "n", lwd = 2  ,xlim=c(0,1.2*max(t_vec)),  ylim = c(0,1.3), main=paste0("Median prevalence estimate with a hunting rate of ", hunt_rates[w]), xlab='time', ylab='observed:true prevalence')
  
  axis( 1, at = t_vec )
  abline(a=1, b=0, col='black', lwd=2, lty='dashed', xlim = c(0,10))
  legend('topright', title = 'Herd Size', box.lty = 0, inset = c(.03,.01), legend = herdSizes, fill = colors)
  for(i in 2:length(ratesToPlot)){
    
    
    lines(x= t_vec, y=replicate(length(t_vec), median(dat[,3*i])/median(dat[,3*i-1])), type='l', col=colors[i], lwd = 2  )
    
  }
  
}
