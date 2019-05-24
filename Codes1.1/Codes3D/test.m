MaxwellDriver3D;
Ez1 = Ez_time;
MaxwellDriver3DLawson;
Ez2 = Ez_time;

figure;
plot(Ez1)
hold on
plot(Ez2)
legend("standard", "lawson")

