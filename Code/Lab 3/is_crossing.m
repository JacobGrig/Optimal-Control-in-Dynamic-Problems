function answer = is_crossing(x1, x2, x3, x4, y1, y2, y3, y4)

gamma = (x1*y3 - x3*y1 - x1*y4 + x4*y1 + x3*y4 - x4*y3)/ ...
    (x1*y3 - x3*y1 - x1*y4 - x2*y3 + x3*y2 + x4*y1 + x2*y4 - x4*y2);

x_cur = x1 + gamma * (x2 - x1);
y_cur = y1 + gamma * (y2 - y1);

answer = (x1 <= x_cur && x_cur <= x2 && x3 <= x_cur && x_cur <= x4 && ...
    y1 <= y_cur && y_cur <= y2 && y3 <= y_cur && y_cur <= y4);

end