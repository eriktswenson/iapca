
*** Swenson - "Building Asymmetry into PCA and Application to ENSO"
*** GrADS script that generates Figure 1 given spatial
*** patterns from sst.DJF.IAPCA.nc and olr.DJF.IAPCA.nc
*** > grads -bpc Fig1.gs

# Author: Erik Swenson (latest revision Sep 2023)

'reinit'

clevs1 = '-1.8 -1.6 -1.4 -1.2 -1 -0.8 -0.6 -0.4 -0.2 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8'
ccols1 = '49 48 47 46 45 44 43 42 41 0 21 22 23 24 25 26 27 28 29'
clevs2 = '5 10 15 20 25 30 35 40 45'

'sdfopen sst.DJF.IAPCA.nc'
'set display color white'; 'c'
'set grads off'; 'set grid off'
'set lon 100 280'; 'set lat -20 20'

'set rgb 21 255 250 170'; 'set rgb 41 200 255 255'
'set rgb 22 255 232 120'; 'set rgb 42 175 240 255'
'set rgb 23 255 192  60'; 'set rgb 43 130 210 255'
'set rgb 24 255 160   0'; 'set rgb 44  95 190 250'
'set rgb 25 255  96   0'; 'set rgb 45  75 180 240'
'set rgb 26 255  50   0'; 'set rgb 46  60 170 230'
'set rgb 27 225  20   0'; 'set rgb 47  40 150 210'
'set rgb 28 192   0   0'; 'set rgb 48  30 140 200'
'set rgb 29 165   0   0'; 'set rgb 49  20 130 190'
'set rgb 56  72  60 200'; 'set rgb 37  55 210  60'

k = 1
while (k<=5)
f = 1; if (k=4); f = -1; endif
'set t 'k
'set parea 1.5 8.5 '11-2*k' '12-2*k'.7'
'set gxout shaded'
'set ylabs 20S||10S||EQ||10N||20N'
'set clevs 'clevs1
'set ccols 'ccols1
'd sst*'f
k = k + 1
endwhile

rc = colorbar(0.7,0,5,0.4)
'set strsiz 0.2'; 'draw string 2.25 0.45 `a0`nC'
'close 1'

'sdfopen olr.DJF.IAPCA.nc'
'set lon 100 280'; 'set lat -20 20'

k = 1
while (k<=5)
f = 1; if (k=4); f = -1; endif
'set t 'k
'set parea 1.5 8.5 '11-2*k' '12-2*k'.7'
'set gxout contour'
'set cthick 4'; 'set clab off'
'set ccolor 37'; 'set clevs 'clevs2
'd olr*-1*'f
'set cthick 4'; 'set clab off'
'set ccolor 56'; 'set clevs 'clevs2
'd olr*'f
k = k + 1
endwhile

'set string 1 c 4'; 'set strsiz 0.18'
'draw string 0.75 10.05 (a)'; 'draw string 0.75 9.65 EOF-1'
'draw string 0.75 8.05 (b)'; 'draw string 0.75 7.65 EOF-2'
'draw string 0.75 6.25 (c)'; 'draw string 0.75 5.85 AEOF-1'; 'draw string 0.75 5.45 positive'
'draw string 0.75 4.25 (d)'; 'draw string 0.75 3.85 AEOF-1'; 'draw string 0.75 3.45 negative'
'draw string 0.75 2.45 (e)'; 'draw string 0.75 2.05 AEOF-1'
'draw string 0.75 1.65 mean'; 'draw string 0.75 1.25 residual'

'gxprint Fig1.pdf'
'quit'


function colorbar (sf,vert,xmid,ymid)

if(sf='');sf=1.0;endif

*  Check shading information
  'query shades'
  shdinfo = result
  if (subwrd(shdinfo,1)='None')
    say "Cannot plot color bar: No shading information"
    return
  endif

*  Get plot size info
  'query gxinfo'
  rec2 = sublin(result,2)
  rec3 = sublin(result,3)
  rec4 = sublin(result,4)
  xsiz = subwrd(rec2,4)
  ysiz = subwrd(rec2,6)
  ylo = subwrd(rec4,4)
  xhi = subwrd(rec3,6)
  xd = xsiz - xhi

  ylolim=0.6*sf
  xdlim1=1.0*sf
  xdlim2=1.5*sf
  barsf=0.8*sf
  yoffset=0.2*sf
  stroff=0.1*sf
  strxsiz=0.22*sf
  strysiz=0.25*sf

*  Decide if horizontal or vertical color bar
*  and set up constants.
  cnum = subwrd(shdinfo,5)
*       logic for setting the bar orientation with user overides
  if (ylo<ylolim | xd>xdlim1)
    vchk = 1
    if(vert = 0) ; vchk = 0 ; endif
  else
    vchk = 0
    if(vert = 1) ; vchk = 1 ; endif
  endif

*       vertical bar
  if (vchk = 1 )

    if(xmid = '') ; xmid = xhi+xd/2 ; endif
    xwid = 0.2*sf
    ywid = 0.5*sf

    xl = xmid-xwid/2
    xr = xl + xwid
    if (ywid*cnum > ysiz*barsf)
      ywid = ysiz*barsf/cnum
    endif
    if(ymid = '') ; ymid = ysiz/2 ; endif
    yb = ymid - ywid*cnum/2
    'set string 1 l 4'
    vert = 1

  else
*       horizontal bar
    ywid = 0.4
    xwid = 0.8

    if(ymid = '') ; ymid = ylo/2-ywid/2 ; endif
    yt = ymid + yoffset
    yb = ymid
    if(xmid = '') ; xmid = xsiz/2 ; endif
    if (xwid*cnum > xsiz*barsf)
      xwid = xsiz*barsf/cnum
    endif
    xl = xmid - xwid*cnum/2
    'set string 1 tc 4'
    vert = 0
  endif

*  Plot colorbar
  'set strsiz 'strxsiz' 'strysiz
  num = 0
  while (num<cnum)
    rec = sublin(shdinfo,num+2)
    col = subwrd(rec,1)
    hi = subwrd(rec,3)

* Note:  Only take NDOT values after decimal point  (LT)
         ndot    = 2
         numdot  = 0
         numdot1 = 1
         while ( numdot1<20 )
         dot = substr(hi,numdot1,1)
         if( dot='.' )
         numdot = numdot1
         endif
         numdot1 = numdot1 + 1
         endwhile
         if( numdot=2 )
             numdot = numdot+2
             hi = substr(hi,1,numdot)
         else
             if( numdot=3 )
                 numdot = numdot+1
                 hi = substr(hi,1,numdot)
             else
                 if( numdot!=0 )
                     hi = substr(hi,1,numdot)
                 endif
             endif
         endif
*        if( numdot!=0 )
*        numdot = numdot+ndot
*        hi = substr(hi,1,numdot)
*        endif

    if (vert)
      yt = yb + ywid
    else
      xr = xl + xwid
    endif

    if(num!=0 & num!= cnum-1)
    'set line 1 1 10'
    'draw rec 'xl' 'yb' 'xr' 'yt
    'set line 'col
    'draw recf 'xl' 'yb' 'xr' 'yt
    even = 0
    if (hi>0); if (math_fmod(num,2)=0); even = 1; endif; endif
    if (hi<0); if (math_fmod(num+1,2)=0); even = 1; endif; endif
    if (num<cnum-1)
      if (vert)
        xp=xr+stroff
        if (even=1); 'draw string 'xp' 'yt' 'hi; endif
      else
        yp=yb-stroff
        if (even=1); 'draw string 'xr' 'yp' 'hi; endif
      endif
    endif
    endif

    if(num = 0 )

      if(vert = 1)

        xm=(xl+xr)*0.5
        'set line 1 1 10'
        'draw line 'xl' 'yt' 'xm' 'yb
        'draw line 'xm' 'yb' 'xr' 'yt
        'draw line 'xr' 'yt' 'xl' 'yt

        'set line 'col
        'draw polyf 'xl' 'yt' 'xm' 'yb' 'xr' 'yt' 'xl' 'yt

      else

        ym=(yb+yt)*0.5
        'set line 1 1 10'
        'draw line 'xl' 'ym' 'xr' 'yb
        'draw line 'xr' 'yb' 'xr' 'yt
        'draw line 'xr' 'yt' 'xl' 'ym

        'set line 'col
       'draw polyf 'xl' 'ym' 'xr' 'yb' 'xr' 'yt' 'xl' 'ym

      endif

    endif

    even = 0
    if (hi>0); if (math_fmod(num,2)=0); even = 1; endif; endif
    if (hi<0); if (math_fmod(num+1,2)=0); even = 1; endif; endif
    if (num<cnum-1)
      if (vert)
         xp=xr+stroff
        if (even=1); 'draw string 'xp' 'yt' 'hi; endif
      else
         yp=yb-stroff
        if (even=1); 'draw string 'xr' 'yp' 'hi; endif
      endif
    endif

    if(num = cnum-1 )

      if( vert = 1)
        'set line 1 1 10'
        'draw line 'xl' 'yb' 'xm' 'yt
        'draw line 'xm' 'yt' 'xr' 'yb
        'draw line 'xr' 'yb' 'xl' 'yb

        'set line 'col
        'draw polyf 'xl' 'yb' 'xm' 'yt' 'xr' 'yb' 'xl' 'yb
      else

        'set line 1 1 10'
        'draw line 'xr' 'ym' 'xl' 'yb
        'draw line 'xl' 'yb' 'xl' 'yt
        'draw line 'xl' 'yt' 'xr' 'ym

        'set line 'col
        'draw polyf 'xr' 'ym' 'xl' 'yb' 'xl' 'yt' 'xr' 'ym

      endif

    endif

    even = 0
    if (hi>0); if (math_fmod(num,2)=0); even = 1; endif; endif
    if (hi<0); if (math_fmod(num+1,2)=0); even = 1; endif; endif
    if (num<cnum-1)
      if (vert)
        xp=xr+stroff
        if (even=1); 'draw string 'xp' 'yt' 'hi; endif
      else
        yp=yb-stroff
       if (even=1); 'draw string 'xr' 'yp' 'hi; endif
      endif
    endif

    num = num + 1
    if (vert); yb = yt;
    else; xl = xr; endif;
  endwhile
return


