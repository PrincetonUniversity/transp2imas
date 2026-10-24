#!/usr/bin/python
"""Regular-expression builders for TRANSP variable and namelist key names with '?' placeholders."""


def make_re_pattern(name:str, nd1:int=0, nd2:int=0, nch1:int=0,
      nch2:int=0, nd:int=0, nch:int=0, dig_last:bool=False, step:int=-1) -> str:
    r'''
      nd1=digit_begin nd2=digit_end. same for char:nch1,nch2.
      nd, nch for fixed number; override nd1,nd2,nch1,nch2
      Used for specified ranges of # of digits and chars
      step=-1 :  add ^ + pattern + $;
      step!=-1 string including 1st '?' after replacing ? by pattern.
      eg: get_var_names_re('VV_PF?', cut='VV_', nd=1, nch1=1, nch2=2).

      #step call
      p = make_re_pattern('NAME\(?,?\)', nd1=1, nd2=3, step=1) => NAME\(\d{1,3}[A-z]{0,0}
      #second of the step call
      make_re_pattern(fr'{p},?\)', nd1=1, nd2=3)  #p is pat from 1st call= \d{1,3}[A-z]{0,0}
                    => NAME\(\d{1,3}[A-z]{0,0},\d{1,3}[A-z]{0,0}\)
      NOTE: escape ().  Error: NAME(?,?)
      For non-regex don't do escape
    '''
    ind = name.find('?') #only 1st position of '?''
    if ind < 0:
        raise Exception(f'Error in pattern build of {name} : missing "?"')

    nd1 = nd if nd>0 else nd1
    nd2 = nd if nd>0 else nd2
    nch1 = nch if nch>0 else nch1
    nch2 = nch if nch>0 else nch2
    pat = f'{name[:ind]}'
    if dig_last:
        pat =  pat + r'[A-z]{' + f'{str(nch1)},{str(nch2)}'+ r'}' + \
             r'\d{' + f'{str(nd1)},{str(nd2)}' + r'}'
    else:
        pat =  pat + r'\d{' + f'{str(nd1)},{str(nd2)}' + r'}' + \
            r'[A-z]{' + f'{str(nch1)},{str(nch2)}' + r'}'

    if step == -1:
        pat = pat + rf'{name[ind+1:]}' + '$'

    return pat


def make_re_pattern_quest(name:str) -> str:
    r'''
     Requires at least one '?' in name
     All '?' are replaced by \w+
     BSTH_?4_P? converted to BSTH_\w+4_P\w+
    '''
    inds = [ind for ind, ch in enumerate(name.lower()) if '?' == ch]
    if len(inds) < 1:
        raise Exception(f'Error in pattern build of {name} : missing "?"')
    pat = rf'{name[:inds[0]]}\w+'
    for i in range(1,len(inds)):
        pat += rf'{name[inds[i-1]+1:inds[i]]}\w+'
    return pat + name[inds[-1]+1:] + '$'


def make_namelist_re_pattern(name:str, nd1:int=0, nd2:int=0, nch1:int=0,
      nch2:int=0, nd:int=0, nch:int=0) -> str:
    #only upto ?,? handled (two occurences)
    nn = name.count('?')
    st = 1 if nn>1 else -1
    pat = make_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch, step=st)
    ii = name.find('?')

    if nn > 1:
        pat = make_re_pattern(fr'{pat}{name[ii+1:]}', nd1, nd2, nch1, nch2, nd, nch)
    return pat
