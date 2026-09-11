# Снимок окна программы для проверки картинки при разработке.
# В проект не входит: это инструмент, а не часть программы.
#
# Запускает molview.exe, при желании щёлкает по точкам клиентской области
# и сохраняет снимок. Клики шлются сообщениями, мышь пользователя не трогается,
# поэтому PostMessage, а не SendMessage: иначе модальное окно (например,
# подтверждение снимка) заблокировало бы сценарий.
#
# Перед запуском нужен путь к DLL библиотеки VTK:
#   $env:PATH = "C:\msys64\mingw64\bin;$env:PATH"
#
# Пример — выбрать пероксид водорода, выключить вращение, сбросить вид:
#   .\tools\snap.ps1 -Exe .\build\molview.exe -Shot out.png -Clicks "60,131;620,745;713,745"

param([string]$Exe, [string]$Shot, [string]$Clicks = "", [int]$Settle = 1500)

Add-Type -AssemblyName System.Drawing
$src = @'
using System;
using System.Runtime.InteropServices;
public class Snp {
  [DllImport("user32.dll")] public static extern bool GetClientRect(IntPtr h, out RECT r);
  [DllImport("user32.dll")] public static extern bool ClientToScreen(IntPtr h, ref POINT p);
  [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr h);
  [DllImport("user32.dll")] public static extern bool SetWindowPos(IntPtr h, IntPtr a, int x, int y, int cx, int cy, uint f);
  [DllImport("user32.dll")] public static extern bool PostMessage(IntPtr h, uint m, IntPtr w, IntPtr l);
  [StructLayout(LayoutKind.Sequential)] public struct RECT { public int Left, Top, Right, Bottom; }
  [StructLayout(LayoutKind.Sequential)] public struct POINT { public int X, Y; }
}
'@
if (-not ("Snp" -as [type])) { Add-Type -TypeDefinition $src }

$p = Start-Process -FilePath $Exe -PassThru
for ($i = 0; $i -lt 40; $i++) {
  Start-Sleep -Milliseconds 250
  $p.Refresh()
  if ($p.HasExited) { "программа завершилась сама, код $($p.ExitCode)"; exit 1 }
  if ($p.MainWindowHandle -ne 0) { break }
}
if ($p.MainWindowHandle -eq 0) { "окно так и не появилось"; $p.Kill(); exit 1 }
$h = $p.MainWindowHandle

Start-Sleep -Milliseconds $Settle
$null = [Snp]::SetWindowPos($h, [IntPtr](-1), 0, 0, 0, 0, 0x0003)
$null = [Snp]::SetForegroundWindow($h)
Start-Sleep -Milliseconds 500

if ($Clicks -ne "") {
  foreach ($pair in $Clicks.Split(";")) {
    $xy = $pair.Split(",")
    $lp = [IntPtr](([int]$xy[1] -shl 16) -bor [int]$xy[0])
    $null = [Snp]::PostMessage($h, 0x0200, [IntPtr]0, $lp)          # WM_MOUSEMOVE
    $null = [Snp]::PostMessage($h, 0x0201, [IntPtr]1, $lp)          # WM_LBUTTONDOWN
    $null = [Snp]::PostMessage($h, 0x0202, [IntPtr]0, $lp)          # WM_LBUTTONUP
    Start-Sleep -Milliseconds 400
  }
  Start-Sleep -Milliseconds 800
}

$p.Refresh()
"ЗАГОЛОВОК: $($p.MainWindowTitle)"
$r = New-Object Snp+RECT; $null = [Snp]::GetClientRect($h, [ref]$r)
$o = New-Object Snp+POINT; $null = [Snp]::ClientToScreen($h, [ref]$o)
$bmp = New-Object System.Drawing.Bitmap $r.Right, $r.Bottom
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.CopyFromScreen($o.X, $o.Y, 0, 0, $bmp.Size)
$bmp.Save($Shot, [System.Drawing.Imaging.ImageFormat]::Png)
"снимок: $Shot"
Start-Sleep -Milliseconds 200
$p.Kill()
