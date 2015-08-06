'''
Miscellaneous functions.
Created on Jul 4, 2014
@author: HSH
'''

def randompartition(lst,n=1):
    '''
    Random partition input list into n equal sub lists.
    '''
    import random
    
    if not isinstance(lst,list):
        raise Exception('First input must be list type.')
    random.shuffle(lst)
    q,r = divmod(len(lst),n)
    indices = [q*i+min(i,r) for i in xrange(n+1)]
    return [lst[indices[i]:indices[i+1]] for i in xrange(n)]
    


def finishNotice(mobile='7186122084'):
    '''Notice via text message when a time consuming script finishes running'''
    import smtplib
    server = smtplib.SMTP('smtp.gmail.com',587)
    server.starttls()
    from_address = 'centralpark2006@gmail.com'
    to_address = mobile + '@txt.att.net'
    server.login(from_address, '11802549')
    server.sendmail(from_address, to_address, 'Python script finished!')
    
    
def errorNotice(mobile='7186122084'):
    '''Notice via text message when a time consuming script raise error'''
    import smtplib
    server = smtplib.SMTP('smtp.gmail.com',587)
    server.starttls()
    from_address = 'centralpark2006@gmail.com'
    to_address = mobile + '@txt.att.net'
    server.login(from_address, '11802549')
    server.sendmail(from_address, to_address, 'Python script error!')
    
def gettestfile(fname,ofname=None,nLines=100,rand=False):
    '''
    Get shorter version of the file for test purposes.
    '''
    import os,random
    
    if not ofname:
        name,ext = os.path.splitext(fname)
        ofname = os.path.join(name+'_text',ext)
    f = open(fname)
    ofile = open(ofname,'w')
    if rand:
        ofile.write(f.readline())
        lines = f.readlines()
        random.shuffle(lines)
        for i in range(nLines):
            ofile.write(lines[i])
    else:
        ofile.write(f.readline())
        for i in range(nLines):
            ofile.write(f.readline())
    f.close()
    ofile.close()
    
    
if __name__ == '__main__':
    gettestfile('/Users/HSH/Roche/Data/all8.txt','/Users/HSH/Roche/20140811/all8_test.txt',rand=True)