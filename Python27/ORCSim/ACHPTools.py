'''
This file holds some functions that don't have any obvious other home
'''
class struct(): pass #an empty class to mimic struct behavior

def dict_to_file(fileobject, dictionary):
    for k,v in dictionary.items():
        fileobject.write('%s = %s\n'%(k,v))
    #end for
#end dict_to_file


def Write2CSV(Class,filename,append=False):
    """
    This function takes in a class and a file pointer
    """
    def OL2Strings(OL):
        head=str(OL[0][0])
        units=str(OL[0][1])
        vals=str(OL[0][2])
        for i in range(1,len(OL)):
            head+=','+str(OL[i][0])
            units+=','+str(OL[i][1])
            vals+=','+str(OL[i][2])
        return head,units,vals
    
    def BuildComponentList(ShapeString,string):
        OutString=string
        for i in range(len(ShapeString.split(','))-1):
            OutString+=','+string
        return OutString
        
#    from Cycle import SecondaryCycleClassORCClass
    from Cycle import ORCClass
    
    if append==False:
        write_or_append = 'w'
    elif append==True:
        write_or_append = 'a'
    
    #the first thing we do is check if the output list is empty. If so, just write a blank line to the file and return
    if Class.OutputList()==[]: #then the list could not be generated because not all required attributes were defined
        if type(filename)!=type('some string'):
            #A file object was passed in, use it
            fP=filename
        else:
            fP=open(filename, write_or_append)
            fP.write('\n')
            fP.close()
        return
    
    # Check if it is an instance of one of the cycle classes - more work required
    # to collect all the component outputs
    # if isinstance(Class,(SecondaryCycleClass,ORCClass)):
    if isinstance(Class,(ORCClass)):
        #Pull the cycle outputs
        head,units,vals=OL2Strings(Class.OutputList())
        headList=[head]
        unitsList=[units]
        valsList=[vals]
        componentList=[BuildComponentList(units,'Cycle')]

        #Loop over the other things that are there
        for item in dir(Class):
            #If the item has an outputList, collect it
            if hasattr(getattr(Class,item),'OutputList'):
                head,units,vals=OL2Strings(getattr(Class,item).OutputList())
                componentList+=[BuildComponentList(units,item)]
                headList+=[head]
                unitsList+=[units]
                valsList+=[vals]
        components=','.join(componentList)
        head=','.join(headList)
        units=','.join(unitsList)
        vals=','.join(valsList)
        IsCycle=True
    else:
        head,units,vals=OL2Strings(Class.OutputList())
        IsCycle=False
    if type(filename)!=type('some string'):
        #A file object was passed in, use it
        fP=filename
    else:
        fP=open(filename, write_or_append)
    
    if append==True:
        fP.write(vals+'\n')
    else:
        if IsCycle==True:
            fP.write(components+'\n')
        fP.write(head+'\n')
        fP.write(units+'\n')
        fP.write(vals+'\n')
    fP.close()








def ValidateFields(d,reqFields,optFields=None):
        '''
        The function checkFields takes in inputs of:
        
        =========   =============================================================
        Variable    Type & Description
        =========   =============================================================
        d           dict of values that are part of structure
        reqFields   list of tuples in the form (fieldname,typepointer,min,max)
        optFields   list of other fieldnames
        =========   =============================================================
        
        required parameters are checked that they 
        * exist
        * can be cast using the typepointer function pointer
        * is within the range (min,max)
        
        if a parameter is on the optional parameters list, it is ok-ed, but not value checked
        
        Additional parameters raise AttributeError
        
        '''
        #make a copy of d
        d=dict(d)
        
        #Required parameters
        for field,typepointer,min,max in reqFields:
            if field in d:
                #See if you can do a type cast using the conversion function pointer
                # You should get the same value back
                assert typepointer(d[field])==d[field],field+': failed type conversion, should be '+str(typepointer)
                #check the bounds if numeric input
                if typepointer in (float,int):
                    assert d[field]>=min and d[field]<=max,field+' (value: %g) not in the range [%g,%g]'%(d[field],min,max)
                #remove field from dictionary of terms left to check if no errors
                del d[field]
            else:
                raise AttributeError('Required field '+field+' not included')
        #Optional parameters (not strictly checked, just checked their existence)
        if optFields!=None:
            for field in optFields:
                if field in d:
                    del d[field]
        assert len(d)==0,'Unmatched fields found: '+str(d.keys())
        
def convert(SourceUnit, TargetUnit):
    if SourceUnit=='in':
        if TargetUnit=='m':
            conversion = 0.0254
        else:
            raise ValueError('Target Unit not supported!')
    elif SourceUnit=='f':
        if TargetUnit=='m':
            conversion = 0.3048
        else:
            raise ValueError('Target Unit not supported!')
    elif SourceUnit=='kg':
        if TargetUnit=='lbm':
            conversion = 2.204622600
        else:
            raise ValueError('Target Unit not supported!')
    elif SourceUnit=='lbm':
        if TargetUnit=='kg':
            conversion = 1/2.204622600
        else:
            raise ValueError('Target Unit not supported!')    
    else:
        raise ValueError('Source Unit not supported!')
    
    return conversion


#The following function reads a file and creates a data structure that holds all the information 
def DataHolderFromFile(filename):
    # class DataHolder: pass #I was going to use a class for data storage, but I cannot easily iterate over property    attributes of an object without iterating over all attributes
    #
    # store = DataHolder() #create an instance of this simple data structure class
    store = {} #use a dictionary to hold all the values in the file by header name
    lines = open(filename).readlines()
    #remove all the newline characters
    i = 0
    for line in lines: #What is the right way to do this? If I modify line, it does not affect the elements of lines.
        lines[i] = lines[i].rstrip('\n')
        i += 1
    
    header = lines[0].split(',') #this gives a list of header names that I will use to define attributes of store
    #create a key in store for each item in header and set the value as an empty list
    for item in header:
        # setattr(store, item, []) #if I used a class, I would do this to create an attribute for each item as an empty list
        store[item] = []
    lines.pop(0) #we no longer need the first row of lines so we remove it
    for line in lines: #loop through all the lines in the file
        line = line.split(',') #split each line into its elements
        
        i = 0
        for item in header:
            # getattr(store, item).append(line[i]) #if I used a class, I could append items to the list by getting the attribute like this
            if i>1:
                store[item].append(float(line[i])) #add an element to the dictionary under the correct key (line has the same number of elements in header
            else:
                store[item].append(line[i])
            i += 1

    #now modify the fluid values in store to make the mixture entries compatible with Refpropp_mix
    i = 0
    for entry in store['fluid']:
        store['fluid'][i] = entry.replace(':', ',') #I had to use a : when I made the file so the fluid entry would not split across columns due to the comma
        i += 1
    return store

def ImportNumericDataColumns(filename):
    store = {} #use a dictionary to hold all the values in the file by header name
    lines = open(filename).readlines()
    #remove all the newline characters
    i = 0
    
    for line in lines: #What is the right way to do this? If I modify line, it does not affect the elements of lines.
        lines[i] = lines[i].rstrip('\n')
        i += 1

    header = lines[0].split(',') #this gives a list of header names that I will use to define attributes of store

    #create a key in store for each item in header and set the value as an empty list
    for item in header:
        # setattr(store, item, []) #if I used a class, I would do this to create an attribute for each item as an empty list
        store[item] = []
    
    lines.pop(0) #we no longer need the first row of lines so we remove it
    for line in lines: #loop through all the lines in the file
        line = line.split(',') #split each line into its elements
        i = 0
    
        for item in header:
            # getattr(store, item).append(line[i]) #if I used a class, I could append items to the list by getting the attribute like this
            store[item].append(float(line[i])) #add an element to the dictionary under the correct key (line has the same number of elements in header
            i += 1

    return store

if __name__=='__main__':
    store_ZRC = DataHolderFromFile('store.csv')