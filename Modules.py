from sympy import *
from math import comb
import numpy as np
import itertools
import os
q=symbols('q')
# '/Users/tinglu/Desktop/Computational_code'

############# data:
A1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
A2 = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
A3 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
A4 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
A5 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
A6 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
A=[A1,A2,A3,A4,A5,A6]

############### dictionaries:
dict2_char2index={
    'a':(0,0),
    'b':(0,1),
    'c':(0,2),
    'd':(1,0),
    'e':(1,1),
    'f':(1,2),
    'g':(2,0),
    'h':(2,1),
    'k':(2,2)
}
dict2_index2char={
    (0,0):'a',
    (0,1):'b',
    (0,2):'c',
    (1,0):'d',
    (1,1):'e',
    (1,2):'f',
    (2,0):'g',
    (2,1):'h',
    (2,2):'k'
}
dict3={
        'a':np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]]),
        'b':np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0]]),
        'c':np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]]),
        'd':np.array([[0, 0, 0], [1, 0, 0], [0, 0, 0]]),
        'e':np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]]),
        'f':np.array([[0, 0, 0], [0, 0, 1], [0, 0, 0]]),
        'g':np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]]),
        'h':np.array([[0, 0, 0], [0, 0, 0], [0, 1, 0]]),
        'k':np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]]),
        'aek':np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        'afh':np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]]),
        'bdk':np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
        'bfg':np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
        'cdh':np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
        'ceg':np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        }
dict_comul={
    'a': (('a','a'),('b','d'),('c','g')),
    'b': (('a','b'),('b','e'),('c','h')),
    'c': (('a','c'),('b','f'),('c','k')),
    'd': (('d','a'),('e','d'),('f','g')),
    'e': (('d','b'),('e','e'),('f','h')),
    'f': (('d','c'),('e','f'),('f','k')),
    'g': (('g','a'),('h','d'),('k','g')),
    'h': (('g','b'),('h','e'),('k','h')),
    'k': (('g','c'),('h','f'),('k','k'))
}


############## methods:
def diff(seq):
    return [seq[i+1] - seq[i] for i in range(len(seq)-1)]

def generator(n, l):
    for combination in itertools.combinations_with_replacement(range(n+1), l-1):
        yield [combination[0]] + diff(combination) + [n-combination[-1]]

def GenCouM(m):
    path=('/Users/tinglu/Desktop/Computational_code/Counting_Matrices/Counting_Matrices_Of_Order_'+str(m))
    if os.path.exists(path)==False:
        os.mkdir(path)
        for a in generator(m,6):
            path1=('/Users/tinglu/Desktop/Computational_code/Counting_Matrices/Counting_Matrices_Of_Order_'+str(m)+'/'+str(a)+'.txt')
            f=open(path1,'w')
            B=a[0]*A[0]+a[1]*A[1]+a[2]*A[2]+a[3]*A[3]+a[4]*A[4]+a[5]*A[5]
            f.write(str(B))  
        return 1
    else:
        return 1
    
def Short2Long(m):
    s=''
    for a in m.split('*'):
            b=a.split('^')
            if len(b)==2:
                for i in range(int(b[1])):
                    s=s+b[0].strip('(').strip(')')
            else:
                s=s+b[0].strip('(').strip(')')
    return s

def MonToCM(m,t):
    Mat=np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    if t=='s':
        for a in m.split('*'):
            b=a.split('^')
            if len(b)==2:
                Mat=Mat+int(b[1])*dict3[b[0].strip('(').strip(')')]
            else:
                for generator in b[0]:
                    Mat=Mat+dict3[generator]
    if t=='l':
        b=list(m)
        for generator in b:
            Mat=Mat+dict3[generator]
    return Mat

def CMtoInd(Mat):
    m=Mat[0][0]+Mat[0][1]+Mat[0][2]
    a=np.min(Mat)
    ind=[a,0,0,a,a,0]
    F=np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
    N=Mat-a*F
    if N[0][0]==0:
        ind[0]=ind[0]+0
        ind[1]=ind[1]+0
        ind[2]=ind[2]+N[2][2]
        ind[3]=ind[3]+N[1][2]
        ind[4]=ind[4]+N[2][1]
        ind[5]=ind[5]+N[1][1]
        return ind
    elif N[0][1]==0:
        ind[0]=ind[0]+N[2][2]
        ind[1]=ind[1]+N[1][2]
        ind[2]=ind[2]+0
        ind[3]=ind[3]+0
        ind[4]=ind[4]+N[1][0]
        ind[5]=ind[5]+N[2][0]
        return ind
    elif N[0][2]==0:
        ind[0]=ind[0]+N[1][1]
        ind[1]=ind[1]+N[2][1]
        ind[2]=ind[2]+N[1][0]
        ind[3]=ind[3]+N[2][0]
        ind[4]=ind[4]+0
        ind[5]=ind[5]+0
        return ind
    elif N[1][0]==0:
        ind[0]=ind[0]+N[2][2]
        ind[1]=ind[1]+N[2][1]
        ind[2]=ind[2]+0
        ind[3]=ind[3]+N[0][1]
        ind[4]=ind[4]+0
        ind[5]=ind[5]+N[0][2]
        return ind
    elif N[1][1]==0:
        ind[0]=ind[0]+0
        ind[1]=ind[1]+N[0][0]
        ind[2]=ind[2]+N[2][2]
        ind[3]=ind[3]+N[2][0]
        ind[4]=ind[4]+N[0][2]
        ind[5]=ind[5]+0
        return ind
    elif N[1][2]==0:
        ind[0]=ind[0]+N[0][0]
        ind[1]=ind[1]+0
        ind[2]=ind[2]+N[0][1]
        ind[3]=ind[3]+0
        ind[4]=ind[4]+N[2][1]
        ind[5]=ind[5]+N[2][0]
        return ind
    elif N[2][0]==0:
        ind[0]=ind[0]+N[1][1]
        ind[1]=ind[1]+N[1][2]
        ind[2]=ind[2]+N[0][1]
        ind[3]=ind[3]+0
        ind[4]=ind[4]+N[0][2]
        ind[5]=ind[5]+0
        return ind
    elif N[2][1]==0:
        ind[0]=ind[0]+N[0][0]
        ind[1]=ind[1]+0
        ind[2]=ind[2]+N[1][0]
        ind[3]=ind[3]+N[1][2]
        ind[4]=ind[4]+0
        ind[5]=ind[5]+N[0][2]
        return ind
    elif N[2][2]==0:
        ind[0]=ind[0]+0
        ind[1]=ind[1]+N[0][0]
        ind[2]=ind[2]+0
        ind[3]=ind[3]+N[0][1]
        ind[4]=ind[4]+N[1][0]
        ind[5]=ind[5]+N[1][1]
        return ind
    

def IndtoStMo(ind,t):
    s=''
    if t=='k' and ind==[1,1,0,1,1,0]:
        s='afhaekbfgcdh'
    elif t=='k' and ind!=[1,1,0,1,1,0]:
        if ind[0]!=0:
            for i in range(ind[0]):
                s=s+'aek'
        if ind[1]!=0:
            for i in range(ind[1]):
                s=s+'afh'
        if ind[2]!=0:
            for i in range(ind[2]):
                s=s+'bdk'
        if ind[3]!=0:
            for i in range(ind[3]):
                s=s+'bfg'
        if ind[4]!=0:
            for i in range(ind[4]):
                s=s+'cdh'
        if ind[5]!=0:
            for i in range(ind[5]):
                s=s+'ceg'
    elif t=='s':
        if ind[0]!=0:
            s=s+'(aek)^'+str(ind[0])+'*'
        if ind[1]!=0:
            s=s+'(afh)^'+str(ind[1])+'*'
        if ind[2]!=0:
            s=s+'(bdk)^'+str(ind[2])+'*'
        if ind[3]!=0:
            s=s+'(bfg)^'+str(ind[3])+'*'
        if ind[4]!=0:
            s=s+'(cdh)^'+str(ind[4])+'*'
        if ind[5]!=0:
            s=s+'(ceg)^'+str(ind[5])
        else:
            s=s[:-1]
    elif t=='l':
        if ind[0]!=0:
            for i in range(ind[0]):
                s=s+'aek'
        if ind[1]!=0:
            for i in range(ind[1]):
                s=s+'afh'
        if ind[2]!=0:
            for i in range(ind[2]):
                s=s+'bdk'
        if ind[3]!=0:
            for i in range(ind[3]):
                s=s+'bfg'
        if ind[4]!=0:
            for i in range(ind[4]):
                s=s+'cdh'
        if ind[5]!=0:
            for i in range(ind[5]):
                s=s+'ceg'
    return s

def generate_dict(m,t):
    dict_1={}
    dict_2={}
    l=[]
    if t=='k':
        for a in generator(m, 6):
            l.append(a)
        for i in range(len(l)):
            dict_1[i]=IndtoStMo(l[len(l)-1-i],t)
            dict_2[IndtoStMo(l[len(l)-1-i],t)]=i
    else:
        for a in generator(m, 6):
            l.append(a)
        for i in range(len(l)):
            dict_1[i]=IndtoStMo(l[len(l)-1-i],t)
            dict_2[IndtoStMo(l[len(l)-1-i],t)]=i
    return dict_1, dict_2

def Source(m,path):
    if os.path.exists(path)==False:
        os.mkdir(path)
    a=factor((-q)**(3*m-2)*(q**2-1)**3*(q**4-1)*(1+q**4-q**2-q**(2*m+2))/(q*(q**(2*m)-1)**2*(q**(2*m+2)-1)**2*(q**(2*m+4)-1)))
    b=factor((-q)**(3*m-2)*(q**2-1)**4*(q**4-1)/((q**(2*m)-1)**2*(q**(2*m+2)-1)**2*(q**(2*m+4)-1)))
    c=factor((-q)**(3*m-1)*(q**2-1)**3*(q**4-1)/((q**(2*m)-1)*(q**(2*m+2)-1)**2*(q**(2*m+4)-1)))
    d=factor((-q)**(3*m)*(q**2-1)**2*(q**4-1)/((q**(2*m+2)-1)**2*(q**(2*m+4)-1)))
    f1=open(path+'/[1,0,0,0,0,'+str(m-1)+'].txt','w')
    f1.write(str(a))
    f1.close()
    f2=open(path+'/[0,1,0,0,0,'+str(m-1)+'].txt','w')
    f2.write(str(b))
    f2.close()
    f3=open(path+'/[0,0,1,0,0,'+str(m-1)+'].txt','w')
    f3.write(str(b))
    f3.close()
    f4=open(path+'/[0,0,0,1,1,'+str(m-2)+'].txt','w')
    f4.write(str(b))
    f4.close()
    f5=open(path+'/[0,0,0,1,0,'+str(m-1)+'].txt','w')
    f5.write(str(c))
    f5.close()
    f6=open(path+'/[0,0,0,0,1,'+str(m-1)+'].txt','w')
    f6.write(str(c))
    f6.close()
    f7=open(path+'/[0,0,0,0,0,'+str(m)+'].txt','w')
    f7.write(str(d))
    f7.close()
    
def Recursive_cg(m,path):
    for i in range(2,m+1):
        f=open(path+str([0,0,0,0,i-1,m-i+1]).replace(' ','')+'.txt','r')
        H=simplify(((i*q**(-2*(m-i)-1)-(i-1)*q**(2*(m-i)+5)-q)/(q**2-1))*eval(str(f.read())))
        f.close()
        if i>2:
            for j in range(2,i):
                c=(1/q-q)**(j-2)*comb(i,j)*q**(-2*(m-i+1))+(q-1/q)**(j-2)*q**(2*j)*comb(i-1,j)*q**(2*(m-i+1))
                f=open(path+str([0,0,0,0,i-j,m-i+j]).replace(' ','')+'.txt','r')
                H=simplify(H-c*eval(str(f.read())))
                f.close()
        f=open(path+str([0,0,0,0,0,m]).replace(' ','')+'.txt','r')
        H=simplify(H-(1/q-q)**(i-2)*q**(-2*(m-i+1))*eval(str(f.read())))
        f.close()
        F=factor(H*(1-q**2)**2/(q**2*(q**(m-i+1)-q**(-m+i-1))**2))
        f=open(path+str([0,0,0,0,i,m-i]).replace(' ','')+'.txt','w')
        f.write(str(F))
        f.close()
        f=open(path+str([0,0,0,i,0,m-i]).replace(' ','')+'.txt','w')
        f.write(str(F))
        f.close()
def Recursive_ccgg(m,path):
    for r in range(1,int(m/2)+1):
        for s in range(m-r):
            f=open(path+str([0,0,0,s,r,m-r-s]).replace(' ','')+'.txt','r')
            H=(-q/(q**2-1))*eval(str(f.read()))
            f.close()
            if s>0:
                for i in range(s):
                    a=(1/q-q)**(i-1)*comb(s+1,i+1)*q**(-2*(m-s))+(q-1/q)**(i-1)*q**(2*i-2)*comb(s,i+1)*q**(2*(m-s+2))
                    f=open(path+str([0,0,0,s-i,r,m-r-s+i]).replace(' ','')+'.txt','r')
                    H=simplify(H-a*eval(str(f.read())))
                    f.close()
            f=open(path+str([0,0,0,0,r,m-r]).replace(' ','')+'.txt','r')
            H=simplify(H-(1/q-q)**(s-1)*q**(-2*(m-s))*eval(str(f.read())))
            f.close()
            F=factor(H*(1-q**2)**2/(q**2*(q**(m-s)-q**(-m+s))**2))
            f=open(path+str([0,0,0,s+1,r,m-r-s-1]).replace(' ','')+'.txt','w')
            f.write(str(F))
            f.close()
            f=open(path+str([0,0,0,r,s+1,m-r-s-1]).replace(' ','')+'.txt','w')
            f.write(str(F))
            f.close()
def find_component(n,string):
    h=[]
    h1=[]
    h2=[]
    h3=[]
    s=list(string)
    if s[0]!='a':
        for i in range(n):
            l1=[]
            for j in range(n):
                if j==i:
                    l1.append(s[0])
                else:
                    l1.append('a')
            h1.append(l1)
    else:
        l1=[]
        for j in range(n):
            l1.append('a')
        h1.append(l1)
    if s[1]!='e':
        for i in range(n):
            l2=[]
            for j in range(n):
                if j==i:
                    l2.append(s[1])
                else:
                    l2.append('e')
            h2.append(l2)
    else:
        l2=[]
        for j in range(n):
            l2.append('e')
        h2.append(l2)
    if s[2]!='k':
        for i in range(n):
            l3=[]
            for j in range(n):
                if j==i:
                    l3.append(s[2])
                else:
                    l3.append('k')
            h3.append(l3)
    else:
        l3=[]
        for j in range(n):
            l3.append('k')
        h3.append(l3)
    for item1 in h1:
        for item2 in h2:
            for item3 in h3:
                a=''
                for p in range(n):
                     a+=item1[p]+item2[p]+item3[p]
                h.append(a)
    return h

def find_component_s(n):
    h=[]
    h1=[]
    h2=[]
    h3=[]
    for i in range(n):
        l1=[]
        for j in range(n):
            if j==i:
                l1.append('b')
            else:
                l1.append('a')
        h1.append(l1)
    for j in range(n-1):
        for k in range(j+1,n):
            l21=[]
            for i in range(n):
                if i==j:
                    l21.append('f')
                elif i==k:
                    l21.append('d')
                else:
                    l21.append('e')
            h2.append(l21)
    for j in range(n-1):
        for k in range(j+1,n):
            l21=[]
            for i in range(n):
                if i==j:
                    l21.append('d')
                elif i==k:
                    l21.append('f')
                else:
                    l21.append('e')
            h2.append(l21)
    for i in range(n):
        l3=[]
        for j in range(n):
            if j==i:
                l3.append('h')
            else:
                l3.append('k')
        h3.append(l3)
    for item1 in h1:
        for item2 in h2:
            for item3 in h3:
                a=''
                for p in range(n):
                     a+=item1[p]+item2[p]+item3[p]
                h.append(a)
    return h

def spec_comp(n,m):
    if n<2:
        return -1
    if m=='afh':
        s1='afk'
        s2='aeh'
        for i in range(n-2):
            s1=s1+'aek'
            s2=s2+'aek'
        s1=s1+'aeh'
        s2=s2+'afk'
        return [s1,s2]
    elif m=='bdk':
        s1='bek'
        s2='adk'
        for i in range(n-2):
            s1=s1+'aek'
            s2=s2+'aek'
        s1=s1+'adk'
        s2=s2+'bek'
        return [s1,s2]
    else:
        return -1

def term_find(string1,string2):
    l=list(string1)
    n=len(l)
    h=find_component(int((n+1)/3),string2)
    t=[]
    for item in h:
        tt=''
        ll=list(item)
        for i in range(n):
            for tu in dict_comul[l[i]]:
                if tu[0]==ll[i]:
                    tt+=tu[1]
                    break
        t.append(tt)
    final=[]
    for j in range(len(h)):
        final.append([h[j],t[j]])
    return final

def term_find_s(string1):
    l=list(string1)
    n=len(l)
    h=find_component_s(int((n+1)/3))
    t=[]
    for item in h:
        tt=''
        ll=list(item)
        for i in range(n):
            for tu in dict_comul[l[i]]:
                if tu[0]==ll[i]:
                    tt+=tu[1]
                    break
        t.append(tt)
    final=[]
    for j in range(len(h)):
        final.append([h[j],t[j]])
    return final

def spec_term_find(string1,string2):
    l=list(string1)
    n=len(l)
    h=spec_comp(int((n+1)/3),string2)
    t=[]
    for item in h:
        tt=''
        ll=list(item)
        for i in range(n):
            for tu in dict_comul[l[i]]:
                if tu[0]==ll[i]:
                    tt+=tu[1]
                    break
        t.append(tt)
    final=[]
    for j in range(len(h)):
        final.append([h[j],t[j]])
    return final

def term_reorder(string,an,kn,t):
    if t=='s':
        c=list(Short2Long(string))
        t=list(IndtoStMo(CMtoInd(MonToCM(string,'s')),'l'))
    if t=='l' or t=='k':
        c=list(str(string))
        t=list(IndtoStMo(CMtoInd(MonToCM(string,'l')),t))
    n=len(c)
    g=[1]
    for i in range(n-1):
        mark=1
        while(mark):
            if c[i] !=t[i]:
                for j in range(i,n):
                    if c[j]==t[i]:
                        m=c[j]
                        c[j]=c[j-1]
                        c[j-1]=m
                        if dict2_char2index[c[j]][0]<dict2_char2index[c[j-1]][0] and dict2_char2index[c[j]][1]==dict2_char2index[c[j-1]][1]:
                            g.append(q)
                        elif dict2_char2index[c[j]][0]>dict2_char2index[c[j-1]][0] and dict2_char2index[c[j]][1]==dict2_char2index[c[j-1]][1]:
                            g.append(q**(-1))
                        elif dict2_char2index[c[j]][0]==dict2_char2index[c[j-1]][0] and dict2_char2index[c[j]][1]<dict2_char2index[c[j-1]][1]:
                            g.append(q)
                        elif dict2_char2index[c[j]][0]==dict2_char2index[c[j-1]][0] and dict2_char2index[c[j]][1]>dict2_char2index[c[j-1]][1]:
                            g.append(q**(-1))
                        elif dict2_char2index[c[j]][0]>dict2_char2index[c[j-1]][0] and dict2_char2index[c[j]][1]>dict2_char2index[c[j-1]][1]:
                            b=c[0:n]
                            b[j-1]=dict2_index2char[(dict2_char2index[c[j-1]][0],dict2_char2index[c[j]][1])]
                            b[j]=dict2_index2char[(dict2_char2index[c[j]][0],dict2_char2index[c[j-1]][1])]
                            g.append(-(q-q**(-1)))
                            a=['']
                            for p in range(n):
                                a[0]+=b[p]
                            c+=a
                        elif dict2_char2index[c[j]][0]<dict2_char2index[c[j-1]][0] and dict2_char2index[c[j]][1]<dict2_char2index[c[j-1]][1]:
                            b=c[0:n]
                            b[j-1]=dict2_index2char[(dict2_char2index[c[j-1]][0],dict2_char2index[c[j]][1])]
                            b[j]=dict2_index2char[(dict2_char2index[c[j]][0],dict2_char2index[c[j-1]][1])]
                            g.append(+(q-q**(-1)))
                            a=['']
                            for p in range(n):
                                a[0]+=b[p]
                            c+=a
                        break
            else:
                mark=0
                
    if c[n-2] !=t[n-2]:
        m=c[n-2]
        c[n-2]=c[n-1]
        c[n-1]=m
        if dict2_char2index[c[n-1]][0]<dict2_char2index[c[n-2]][0] and dict2_char2index[c[n-1]][1]==dict2_char2index[c[n-2]][1]:
            c.append(q)
        elif dict2_char2index[c[n-1]][0]>dict2_char2index[c[n-2]][0] and dict2_char2index[c[n-1]][1]==dict2_char2index[c[n-2]][1]:
            c.append(q**(-1))
        elif dict2_char2index[c[n-1]][0]==dict2_char2index[c[n-2]][0] and dict2_char2index[c[n-1]][1]<dict2_char2index[c[n-2]][1]:
            c.append(q)
        elif dict2_char2index[c[n-1]][0]==dict2_char2index[c[n-2]][0] and dict2_char2index[c[n-1]][1]>dict2_char2index[c[n-2]][1]:
            c.append(q**(-1))
        elif dict2_char2index[c[n-1]][0]>dict2_char2index[c[n-2]][0] and dict2_char2index[c[n-1]][1]>dict2_char2index[c[n-2]][1]:
            b=c[0:n]
            b[n-2]=dict2_index2char[(dict2_char2index[c[n-2]][0],dict2_char2index[c[n-1]][1])]
            b[n-1]=dict2_index2char[(dict2_char2index[c[n-1]][0],dict2_char2index[c[n-2]][1])]
            g.append(-(q-q**(-1)))
            a=['']
            for p in range(n):
                a[0]+=b[p]
            c+=a
        elif dict2_char2index[c[n-1]][0]<dict2_char2index[c[n-2]][0] and dict2_char2index[c[n-1]][1]<dict2_char2index[c[n-2]][1]:
            b=c[0:n]
            b[n-2]=dict2_index2char[(dict2_char2index[c[n-2]][0],dict2_char2index[c[n-1]][1])]
            b[n-1]=dict2_index2char[(dict2_char2index[c[n-1]][0],dict2_char2index[c[n-2]][1])]
            g.append(+(q-q**(-1)))
            a=['']
            for p in range(n):
                a[0]+=b[p]
            c+=a
    else:
        f=[]
        f.append(c[0])
        s=len(c)
        if s>n:
            for t in range(1,n):
                f[0]+=c[t]
            for t in range(n,s):
                f.append(c[t])
        else:
            for t in range(1,n):
                f[0]+=c[t]
    ap=[]
    if type(an)==int:
        for u in range(len(f)):
            a_count=0
            aal=list(f[u])
            for k in aal:
                if k=='a':
                    a_count+=1
            if a_count<an:
                ap.append(u)  
    else:
        return f,g
    kp=[]
    if type(kn)==int:
        for u in range(len(f)):
            k_count=0
            kll=list(f[u])
            for k in kll:
                if k=='k':
                    k_count+=1
            if k_count<kn:
                kp.append(u)
    fp=ap
    for item in kp:
        if item in fp:
            continue
        else:
            fp.append(item)
    fp.sort(reverse=True)
    gp=[]
    for j in range(len(g)):
        if len(str(g[j]))>4:
            gp.append(j)
    gp.sort()
    for i in fp:
        f.pop(i)
    for y in range(len(fp)):
        g.pop(gp[fp[y]-1])
    return f, g

def comp_coeff(l1,l2):
    l3=[]
    l4=[]
    l5=[]
    c=1
    n=len(l1)
    for i in range(len(l2)):
        if len(str(l2[i]))>4:
            l3.append(i)
        else:
            l4.append(i)
    for item in l4:
        c=c*l2[item]
    l5.append(c)
    for item in l3:
        c=1
        for i in range(item):
            if len(str(l2[i]))<4:
                c=c*l2[i]
        c=c*l2[item]
        l5.append(c)
    return l1, l5

def term_expand(string,an,kn,t):
    [f,g]=term_reorder(string,an,kn,t)
    [l1,l2]=comp_coeff(f,g)
    a=0
    if t=='s':
        z='l'
    else:
        z=t
    mark=1
    while(mark):
        for i in range(len(l1)):
            if l1[i]!=IndtoStMo(CMtoInd(MonToCM(l1[i],'l')),z):
                [h,n]=term_reorder(l1[i],an,kn,z)
                [l3,l4]=comp_coeff(h,n)
                l1.pop(i)
                l1=l1+l3
                for j in range(len(l4)):
                    l4[j]=l4[j]*l2[i]
                l2.pop(i)
                l2=l2+l4
            else:
                a=a+1
        if a==len(l1):
            mark=0
        else:
            a=0
    return l1,l2

def term_expansion(string,an,kn,t,dict3_int2string,dict3_string2int,p):
    [term,coef]=term_expand(string,an,kn,t)
    h=[0]*len(dict3_string2int)
    for i in range(len(term)):
        h[dict3_string2int[term[i]]]=h[dict3_string2int[term[i]]]+coef[i]
    c=''
    w=[]
    for i in range(len(h)):
        h[i]=simplify(h[i])
        if h[i]!=0:
            w.append(i)
    n=len(w)-1
    for i in range(len(w)):
        if h[w[i]]==1 and len(w)>1:
            c=c+dict3_int2string[w[i]]+'+'
        elif h[w[i]]==1:
            c=c+dict3_int2string[w[i]]
        elif i==n:
            c=c+'('+str(h[w[i]])+')'+'*'+dict3_int2string[w[i]]
        else:
            c=c+'('+str(h[w[i]])+')'+'*'+dict3_int2string[w[i]]+'+'
    if p==1:
        print(string,'=',c)
    return h

def generate_linear_relation(list2,string,order,dict_1,dict_2):
    if string=='afh':
        a=1
        b=0
    if string=='bdk':
        a=0
        b=1
    z=term_find(IndtoStMo(list2,'l'),string)
    t=[0]*len(dict_1)
    for item in z:
        c1=term_expansion(item[0],order-b,order-a,'l',dict_1,dict_2,0)
        c2=term_expansion(item[1],0,0,'l',dict_1,dict_2,0)
        for i in range(len(dict_1)):
            t[i]=factor(t[i]+c1[dict_2[IndtoStMo([order-1,a,b,0,0,0],'l')]]*c2[i])
    t[dict_2[IndtoStMo(list2,'l')]]=t[dict_2[IndtoStMo(list2,'l')]]+order*q
    return t

def Linear_Haar_Compute(list1,list2,string,order,dict_1,dict_2,path):
    t=generate_linear_relation(list2,string,order,dict_1,dict_2)
    r=0
    for i in range(len(dict_1)):
        if i!=dict_2[IndtoStMo(list1,'l')] and t[i]!=0:
            f=open(path+str(CMtoInd(MonToCM(dict_1[i],'l'))).replace(' ','')+'.txt','r')
            H=eval(str(f.read()))
            f.close()
            r=simplify(r+t[i]*H)
    d=factor(-r/t[dict_2[IndtoStMo(list1,'l')]])
    return d
