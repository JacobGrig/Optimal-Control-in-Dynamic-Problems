function [X, Y] = delete_cross(X, Y)

X_new = X;
Y_new = Y;

for index = 3 : length(X)
    
    x1 = X(index);
    y1 = Y(index);
    
    x2 = X(index - 1);
    y2 = Y(index - 1);
    
    for cur_idx = 1 : (index - 1)
        
        x3 = X(cur_idx);
        y3 = Y(cur_idx);
        
        x4 = X(cur_idx + 1);
        y4 = Y(cur_idx + 1);
        
        if (is_crossing(x1, x2, x3, x4, y1, y2, y3, y4)
            
            
            break;
            
        end
        
        
    end
    
end

end