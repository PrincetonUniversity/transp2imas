#!/usr/bin/python3
#
class imas_connection:
    def __init__(self, shot, run, backend, occur, user):
        self.shot = shot
        self.run = run
        # imas_env()
        # self.ids = either create or open existing one
        self.ids=1
        print("we are in imas_connection", self.shot, self.run)
        
class netcdf_connection:
    def __init__(self, filename):
        self.filename = filename
        #self.nc = file handler for netcdf
        self.nc = 2
        print("We are in netcdf_connection")
      
class nc2imas(imas_connection, netcdf_connection):
    def __init__(self, shot, run, backend, occur, user, filename):
        super().__init__(shot, run, backend, occur, user)
        super(imas_connection, self).__init__(filename)
        print("We are in nc2imas")
   
   # def ...fill_values will be here

class imas_summary(nc2imas):
    def __init__(self, ncimas:nc2imas):
        # self.ids = ncimas.ids.summary
        self.nc = ncimas.nc
        print("We are in imas_summary", self.nc)

nc = nc2imas(123124, 1, 1, 0, 'pankin', 'TRANSP.CDF')
sum = imas_summary(nc)


