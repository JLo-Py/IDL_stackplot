function line_perp, x, a, b, perp=perp, point=point, m=m, n=n

if not(keyword_set(point)) and keyword_set(perp) then begin

    point=[0,0]
    print, 'The perpendicular line will pass through [0,0]'
    
endif 
    
if keyword_set(perp) then begin 
    
    m=-float(1)/a
    n=point[1]-m*point[0]     
    
    y=m*x+n
        
endif else begin 
    
    y=a*x+b
    
endelse 
    
return, y

end 
