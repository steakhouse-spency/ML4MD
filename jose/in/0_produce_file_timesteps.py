'''
Created on Nov 1, 2015

@author: Jose_Cobena
'''

def main():

    #name_file_ccf = 'time-steps-ccf.times'
    name_file_isf = 'time-steps-isf.times'
    #name_file_vac = 'time-steps-vac.times'

#---------------------------------------------------------------
    with open(name_file_isf, 'w') as time_steps_file:
        time_steps = []

        for i in range(2,11,2):
            time_steps.append(i)
        
        for i in range(20,101,10):
            time_steps.append(i)
        
        for i in range(200,1001,100):
            time_steps.append(i)
            
        for i in range(2000,10001,1000):
            time_steps.append(i)
        
        for i in range(20000,100001,10000):
            time_steps.append(i)
        
        for i in range(200000,3100001,100000):
            time_steps.append(i)
        
        for j in time_steps:
            print(j, file = time_steps_file)


if __name__ == "__main__": main()
