UAF_USRP
========
  
Software for using the USRP as an ionosonde  

General operation:  
1. Start the ionosonde server
$ UAF_USRP/ionosonde_server/ionosonde_server  
  
2. Start a control program  
$ UAF_USRP/control_programs/swept_freq/swept_freq --[OPTIONS]  
  
3. View the results!  
$ UAF_USRP/utils/h5view.py --[DATE] --[TIME]  
  
See the UAF_USRP/docs directory for more information.

