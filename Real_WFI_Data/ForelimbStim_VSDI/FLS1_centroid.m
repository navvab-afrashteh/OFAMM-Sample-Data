centro = [34,48,56;
35,48,53;
36,51,49;
37,59,45;
38,61,44;
39,58,45;
40,53,46;
41,53,48;
42,51,49;];

points = centro(:,2:3);

dP = diff(points);
dT= diff(centro(:,1));
distances = sqrt(dP(:,1).^2 + dP(:,2).^2);
speeds = distances./dT;

figure(10);clf;
plot(speeds);