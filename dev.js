const { exec, spawn } = require('child_process');
const path = require('path');

const file = process.argv[2];
if (!file) {
    console.error('请提供要编译的 C++ 文件，例如：npm run cpp -- b.cpp');
    process.exit(1);
}

if (path.extname(file) !== '.cpp') {
    console.error('请输入一个 .cpp 文件');
    process.exit(1);
}

const baseName = path.basename(file, '.cpp');
const isWin = process.platform === 'win32';
const output = isWin ? `${baseName}.exe` : `${baseName}.out`;

const compileCmd = `g++ -std=c++23 ${file} -o ${output}`;
console.log(`编译命令：${compileCmd}`);

exec(compileCmd, (compileErr, compileStdout, compileStderr) => {
    if (compileErr) {
        console.error(`编译失败:\n${compileStderr}`);
        return;
    }
    console.log('编译成功！');

    // 使用 spawn 来运行编译后的程序，并允许交互输入
    console.log(`运行命令：${output}`);
    const runProcess = spawn(output, [], { shell: true, stdio: 'inherit' });

    runProcess.on('close', (code) => {
        console.log(`程序退出，退出码 ${code}`);
    });
});
